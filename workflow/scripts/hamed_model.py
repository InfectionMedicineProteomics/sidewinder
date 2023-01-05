#! /usr/bin/env python3


import os
from pathlib import Path

import click
from Bio.PDB import PDBIO, PDBParser

from utils.cheetah import xlvalidation as xlv, megadock_run as mdr
#from .xlgenerator import xlgenerator


def modeling(complex_pdb: Path,
             top_xls: Path,
             output_dir: Path,
             dock_file: Path, cut_off: int, n_models: int, n_filters: int):

    parser = PDBParser()

    structure_input = parser.get_structure('INPS', complex_pdb)

    model = structure_input[0]

    partner_A_pdb = model['A']

    partner_B_pdb = model['B']

    pA_out = output_dir / 'pdb_A.pdb'

    pB_out = output_dir / 'pdb_B.pdb'

    # Rewriting pdb_A and pdb_B to make sure the ending is correct (TER and END)
    # megadock cannot produce correct models with wrong ending
    io = PDBIO()

    io.set_structure(partner_A_pdb)

    io.save(str(pA_out))

    io.set_structure(partner_B_pdb)

    io.save(str(pB_out))

    score_t_list = []

    score_normal_list = []

    structure_list = []

    xls_list = []

    xls_list_out = []

    out_PDBs_list = []

    model_dir = output_dir / "lr_models"

    model_dir.mkdir(parents=True, exist_ok=True)

    job_output = model_dir / 'lr_model'

    # This counter will stop the calculation if after 10 rounds (can be user specified) it can't
    # find a new pose justifying more XLs. Then we report the pose with MAX number of XLs.
    break_counter = 1
    # num_of_models = 10 #for test
    # n_filters = 5 #for test

    #top_xls_name = 'top_model_XLs.txt'
    #top_model_XLs_file = open(top_xls_name, 'w')

    # Global docking (low resolution)
    modeling_counter = 1

    print("Modeling is started ...\n")

    while break_counter <= 100:

        print("Round number ", break_counter)

        dock_structure = mdr.megadock_run(pA_out,
                                          pB_out,
                                          dock_file,
                                          modeling_counter, output_dir)

        modeling_counter += 1

        if dock_structure is not None:

            (score_t,
             score_normal,
             xl_below_cut) = xlv.xlvalidation(dock_structure,
                                              str(top_xls), cut_off)

            # keeping top structures according to the "n_filters"
            if len(score_t_list) > n_filters:
                min_score = min(score_t_list)
                index_min = score_t_list.index(min_score)
                if score_t > min_score:
                    score_t_list[index_min] = score_t
                    score_normal_list[index_min] = score_normal
                    structure_list[index_min] = dock_structure
                    xls_list[index_min] = xl_below_cut

            else:
                score_t_list.append(score_t)
                score_normal_list.append(score_normal)
                structure_list.append(dock_structure)
                xls_list.append(xl_below_cut)

            break_counter += 1

    # storing all the top models
    for num_pos, struct in enumerate(structure_list):
        pdb_name = str(job_output)+"_"+str(num_pos)+".pdb"
        xl_file_name = str(job_output)+"_"+str(num_pos)+"_xls.txt"
        xls_list_out.append(xl_file_name)
        out_PDBs_list.append(pdb_name)

        io = PDBIO()
        io.set_structure(struct)
        io.save(pdb_name)

        with open(xl_file_name, 'w') as f:
            for item in xls_list[num_pos]:
                f.write("%s\n" % item)

        #xlgenerator(struct, "A_B", cut_off, out_name)

    return xls_list_out, out_PDBs_list, score_t_list


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
@click.option('--top_xls',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--dock_file',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--n_models',
              type=int,
              default=10,
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
              output_dir, top_xls, dock_file, n_models, n_filters, cut_off):
    """
    """

    (xl_files_list,
     pdb_files_list,
     score_t_list) = modeling(complex_pdb,
                              top_xls,
                              output_dir,
                              dock_file, cut_off, n_models, n_filters)

    # TODO:
    # Rename the top model and XL files...
    top_model_xl = xl_files_list[score_t_list.index(max(score_t_list))]
    top_model = pdb_files_list[score_t_list.index(max(score_t_list))]

    command1 = "cp %s %s/best_model.pdb" % (top_model, output_dir)
    os.system(command1)

    command2 = "cp %s %s/best_model_xl.txt" % (top_model_xl, output_dir)
    os.system(command2)


if __name__ == "__main__":

    run_model()
