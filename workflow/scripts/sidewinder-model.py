#!/usr/bin/env python3
"""Cheetah-MS docking module.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Save euclidean distances to variable and dump with output


import os
from pathlib import Path
from typing import List, Tuple

import click
from Bio.PDB import PDBIO, PDBParser

from utils import xlvalidation as xlv, megadock_run as mdr


def modeling(pdb_a: Path,
             pdb_b: Path,
             top_xls_file: Path,
             output_dir: Path,
             dock_file: Path,
             cut_off: int,
             n_models: int,
             n_filters: int) -> Tuple[List[List[str]],
                                      List[Path], List[int], List[float]]:
    """
    Performs protein-protein docking using Megadock, validates models with cross-links,
    and filters top models based on the number of satisfied cross-links and score.

    This function conducts the core protein-protein docking and filtering steps.
    It parses the input PDB, performs docking with Megadock, validates the generated
    models with cross-links, and filters the top models based on user-specified criteria.

    Args:
        complex_pdb (Path): Path to the input PDB file containing the protein complex.
        top_xls_file (Path): Path to the file containing top cross-links.
        output_dir (Path): Path to the desired output directory for results.
        dock_file (Path): Path to the docking parameter file for Megadock.
        cut_off (int): Distance cutoff for validating cross-links.
        n_models (int): Number of models to generate in Megadock.
        n_filters (int): Number of models to keep after filtering based on XL hits and score.

    Returns:
        Tuple[List[List[str]], List[Path], List[int], List[float]]:
            - List[List[str]]: List of lists containing validated cross-links (XLs) for each model.
            - List[Path]: List of Path objects to the output PDB files for the top models.
            - List[int]: List containing the number of satisfied XLs for each top model.
            - List[float]: List containing the score for each top model.

    Raises:
        ValueError: If invalid paths are provided for input files or output directory.
        FileNotFoundError: If specified input files cannot be found.
    """
    # parser = PDBParser()

    # structure_input = parser.get_structure('INPS', complex_pdb)

    # model = structure_input[0]

    # partner_A_pdb = model['A']

    # partner_B_pdb = model['B']

    # pA_out = output_dir / 'pdb_A.pdb'

    # pB_out = output_dir / 'pdb_B.pdb'

    top_xls = []

    with open(top_xls_file, 'r') as f:

        for line in f.read().splitlines():

            top_xls.append(line.strip())

    # Rewriting pdb_A and pdb_B to make sure the ending is correct (TER and END)
    # megadock cannot produce correct models with wrong ending.
    # io = PDBIO()

    # io.set_structure(partner_A_pdb)

    # io.save(str(pA_out))

    # io.set_structure(partner_B_pdb)

    # io.save(str(pB_out))

    n_xls_list = []

    score_list = []

    structure_list = []

    xls_list = []

    xls_list_out = []

    out_pdb_list = []

    model_dir = output_dir / 'lr_models'

    model_dir.mkdir(parents=True, exist_ok=True)

    job_output = model_dir / 'lr_model'

    # Find a new pose justifying more XLs. Then we report the pose with MAX number of XLs.
    break_counter = 1

    # Global docking (low resolution)
    modeling_counter = 1

    while break_counter <= n_models:

        dock_structure = mdr.megadock_run(pdb_a,
                                          pdb_b,
                                          dock_file,
                                          modeling_counter, output_dir)

        modeling_counter += 1

        if dock_structure is not None:

            (score_t,
             score_normal,
             xl_below_cut) = xlv.xlvalidation(dock_structure,
                                              top_xls, cut_off)

            # Keeping top structures according to the "n_filters"
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

    return (xls_list_out, out_pdb_list, n_xls_list, score_list)

def get_best_model_idx(xl_hits: List[int], xl_scores: List[float]) -> int:
    """Identifies the index of the model with the most satisfied cross-links (XLs) and highest score.

    This function determines the index of the best model among a set of generated models.
    It considers two criteria: the number of satisfied XLs (XL hits) and the score of the model.

    The function prioritizes models with the highest number of satisfied XLs.
    If multiple models have the same highest number of XL hits, the model with the highest score is chosen.

    Args:
        xl_hits (List[int]): List containing the number of satisfied XLs for each model.
        xl_scores (List[float]): List containing the score for each model.

    Returns:
        int: The index of the model with the most satisfied XLs and highest score.
    """
    most_xls = max(xl_hits)

    if all(n == most_xls for n in xl_hits):

        best_index = xl_scores.index(max(xl_scores))

    else:

        best_index = xl_hits.index(most_xls)

    return best_index


@click.command()
@click.version_option(version='1.0a')
@click.option('--pdb_a',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='')
@click.option('--pdb_b',
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
              default=10,
              help='')
@click.option('--n_filters',
              type=int,
              default=3,
              help='')
@click.option('--cut_off',
              type=int,
              default=32,
              help='')
def run_model(pdb_a: Path,
              pdb_b: Path,
              output_dir: Path,
              top_xls_file: Path,
              dock_file: Path,
              n_models: int = 10,
              n_filters: int = 3, cut_off: int = 32):
    """Runs the protein-protein docking and model filtering process.

    This function serves as the command-line interface for the docking workflow.
    It handles user inputs, calls the necessary functions, and manages output files.

    Args:
        complex_pdb (Path): Path to the input PDB file containing the protein complex.
        output_dir (Path): Path to the desired output directory for results.
        top_xls_file (Path): Path to the file containing top cross-links.
        dock_file (Path): Path to the docking parameter file for Megadock.
        n_models (int, optional): Number of models to generate in Megadock. Defaults to 10.
        n_filters (int, optional): Number of models to keep after filtering. Defaults to 3.
        cut_off (int, optional): Distance cutoff for validating cross-links. Defaults to 32.

    Raises:
        ValueError: If invalid paths are provided for input files or output directory.
        FileNotFoundError: If specified input files cannot be found.
    """
    n_filters -= 1

    (xl_files,
     pdb_files,
     xl_hits,
     xl_scores) = modeling(pdb_a,
                           pdb_b,
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
