#!/usr/bin/env python3

# TODO:
# - Save euclidean distances to variable and dump with output


import os
from pathlib import Path
from typing import List

import click
from Bio.PDB import PDBIO, PDBParser

from utils.cheetah import xlvalidation as xlv, megadock_run as mdr
#from .xlgenerator import xlgenerator


def modeling(complex_pdb: Path,
             top_xls_file: Path,
             output_dir: Path,
             dock_file: Path, cut_off: int, n_models: int, n_filters: int):

    parser = PDBParser()

    structure_input = parser.get_structure('INPS', complex_pdb)

    model = structure_input[0]

    partner_A_pdb = model['A']

    partner_B_pdb = model['B']

    pA_out = output_dir / 'pdb_A.pdb'

    pB_out = output_dir / 'pdb_B.pdb'

    top_xls = []

    with open(top_xls_file, 'r') as f:

        for line in f.read().splitlines():

            top_xls.append(line.strip())

    # Rewriting pdb_A and pdb_B to make sure the ending is correct (TER and END)
    # megadock cannot produce correct models with wrong ending
    io = PDBIO()

    io.set_structure(partner_A_pdb)

    io.save(str(pA_out))

    io.set_structure(partner_B_pdb)

    io.save(str(pB_out))

    n_xls_list = []

    score_list = []

    structure_list = []

    xls_list = []

    xls_list_out = []

    out_pdb_list = []

    model_dir = output_dir / 'lr_models'

    model_dir.mkdir(parents=True, exist_ok=True)

    job_output = model_dir / 'lr_model'

    # find a new pose justifying more XLs. Then we report the pose with MAX number of XLs.
    break_counter = 1

    # Global docking (low resolution)
    modeling_counter = 1

    while break_counter <= n_models:

        dock_structure = mdr.megadock_run(pA_out,
                                          pB_out,
                                          dock_file,
                                          modeling_counter, output_dir)

        modeling_counter += 1

        if dock_structure is not None:

            (score_t,
             score_normal,
             xl_below_cut) = xlv.xlvalidation(dock_structure,
                                              top_xls, cut_off)

            # keeping top structures according to the "n_filters"
            if len(n_xls_list) > n_filters:

                min_score = min(n_xls_list)

                index_min = n_xls_list.index(min_score)

                if score_t > min_score:

                    n_xls_list[index_min] = score_t

                    score_list[index_min] = score_normal

                    structure_list[index_min] = dock_structure

                    xls_list[index_min] = xl_below_cut

            else:

                n_xls_list.append(score_t)

                score_list.append(score_normal)

                structure_list.append(dock_structure)

                xls_list.append(xl_below_cut)

            break_counter += 1

    # Storing all the top models.
    for num_pos, struct in enumerate(structure_list):

        pdb_name = str(job_output)+"_"+str(num_pos)+".pdb"

        xl_file_name = str(job_output)+"_"+str(num_pos)+"_xls.txt"

        xls_list_out.append(xl_file_name)

        out_pdb_list.append(pdb_name)

        io = PDBIO()

        io.set_structure(struct)

        io.save(pdb_name)

        with open(xl_file_name, 'w') as f:

            for item in xls_list[num_pos]:

                f.write(f'{item}\n')

    score_file = model_dir / 'scores.csv'

    with open(score_file, 'w') as f:

        f.write(f'model,xl_hits,score\n')

        for idx, model_path in enumerate(out_pdb_list):

            model = Path(model_path).stem

            f.write(f'{model},{n_xls_list[idx]},{score_list[idx]}\n')

    return xls_list_out, out_pdb_list, n_xls_list, score_list

def get_best_model_idx(xl_hits: List[int], xl_scores: List[float]):

    most_xls = max(xl_hits)

    if all(n == most_xls for n in xl_hits):

        best_index = xl_scores.index(max(xl_scores))

    else:

        best_index = xl_hits.index(most_xls)

    return best_index


@click.command()
@click.version_option(version='1.0a')
@click.option('--complex_pdb',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--output_dir',
              '-o',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--top_xls_file',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--dock_file',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--n_models',
              type=int,
              default=100,
              help='')
@click.option('--n_filters',
              type=int,
              default=2,
              help='')
@click.option('--cut_off',
              type=int,
              default=35,
              help='')
def run_model(complex_pdb,
              output_dir,
              top_xls_file, dock_file, n_models, n_filters, cut_off):
    """
    """

    (xl_files,
     pdb_files,
     xl_hits,
     xl_scores) = modeling(complex_pdb,
                           top_xls_file,
                           output_dir, dock_file, cut_off, n_models, n_filters)

    # TODO:
    # Rename the top model and XL files...
    idx = get_best_model_idx(xl_hits, xl_scores)

    top_model = pdb_files[idx]

    top_model_xls = xl_files[idx]

    command1 = f'cp {top_model} {output_dir}/best_model.pdb'

    os.system(command1)

    command2 = f'cp {top_model_xls} {output_dir}/best_model_xls.txt'

    os.system(command2)


if __name__ == "__main__":

    run_model()
