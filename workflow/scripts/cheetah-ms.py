#!/usr/bin/env python3
"""
"""


from glob import glob
from pathlib import Path

import click
#import pyrosetta
#import pyrosetta.rosetta as rosetta
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta import (
    init,
    Pose,
    pose_from_pdb,
    PyJobDistributor,
    create_score_function
)

from utils.cheetah import (
    rosettaxl,
    rosettaxlv
)
from utils import file_handling


@click.command()
@click.version_option(version='0.2')
@click.option(
    '--input_models_dir',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Input models directory'
)
@click.option(
    '--output_dir',
    '-o',
    required=True,
    type=click.Path(writable=True, path_type=Path),
    help='Output directory'
)
# @click.option(  # This file should be produced by cheetah-ms.py.
#     '--top_xl_file',
#     '-t',
#     required=True,
#     type=click.Path(exists=True, path_type=Path),
#     help='Top ranked cross-links (in Kojak format) based on MS1 features'
# )
@click.option(
    '--cutoff',
    type=int,
    default=35,
    help=(
        f'Cutoff for Euclidean distance measurements (Ånström) between'
        f'cross-linked peptides'
    )
)
@click.option(
    '--partner_chains',
    type=str,
    default='A_B',
    help='Interaction partners chain names in structure file'
)
@click.option(
    '--n_top_models',
    type=int,
    default=3,
    help=''
)
def main(input_models_dir,
         output_dir, top_xl_file, cutoff, partner_chains, n_top_models):
    """...

    ...

    Parameters
    ----------


    Returns
    -------

    """

    init()  # Initiates Rosetta.

    input_models = file_handling.case_ins_ext_glob(

        directory=input_models_dir,
        extensions=('.pdb', '.mmcif')

    )

    score_t_list = []
    score_normal_list = []
    pose_list = []
    out_all_xls_list = []
    out_xls_list = []
    xl_below_cutoff_list = []

    for model_path in input_models:

        dock_pose = pose_from_pdb(model_path)

        score_t, score_normal, xl_below_cutoff = rosettaxlv.rosettaxlv(

            dock_pose,
            top_xl_file,
            cutoff

        )

        if len(score_t_list) > n_top_models:

            min_score = min(score_t_list)

            index_min = score_t_list.index(min_score)

            if score_t > min_score:

                score_t_list[index_min] = score_t

                score_normal_list[index_min] = score_normal

                pose_list[index_min] = dock_pose

                xl_below_cutoff_list[index_min] = xl_below_cutoff

        else:

            score_t_list.append(score_t)

            score_normal_list.append(score_normal)

            pose_list.append(dock_pose)

            xl_below_cutoff_list.append(xl_below_cutoff)

    for num_pos, struct in enumerate(pose_list):

        pdb_name = Path(f'lr_model_{str(num_pos)}')

        full_path = str(output_dir / pdb_name)  # Put in proper file check!

        files = [

            x for y in ('*.fasc', '*.pdb', '*.txt')

            for x in output_dir.glob(f'{pdb_name}{y}')

        ]

        if len(files) > 0:

            for f in files:

                Path(f).unlink()

        out_all_name = output_dir.joinpath(f'{pdb_name}_all_xls.txt')

        out_all_xls_list.append(str(out_all_name))

        scorefxn_low = create_score_function('interchain_cen')

        jd = PyJobDistributor(full_path, 1, scorefxn_low)

        struct.pdb_info().name(full_path + '_fa')

        jd.output_decoy(struct)

        rosettaxl.rosettaxl(struct, partner_chains, cutoff, str(out_all_name))

        out_name = output_dir.joinpath(f'{pdb_name}_xls.txt')

        out_xls_list.append(out_name)

        with open(out_name, 'w') as f:

            XLs = xl_below_cutoff_list[num_pos]

            for XL in XLs:

                print(XL, file=f)


if __name__ == "__main__":

    main()
