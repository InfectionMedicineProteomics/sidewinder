#!/usr/bin/env python3
"""Cheetah-MS docking module.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


import __main__
import csv
import pickle
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import List
__main__.pymol_argv = ['pymol', '-qc']

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol as pm
import seaborn as sns
from abnumber import Chain
from Bio import SeqIO
from Bio.PDB import PDBParser
from matplotlib.patches import Rectangle
from pandas import DataFrame, Series
from scipy.integrate import simps
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import distance as Distance

from utils import dataframe


@dataclass
class DataLoader:
    """Load Sidewinder output database as Pandas DataFrame"""

    df = None  # Final concatenated DataFrame

    def load_data(self,
                  project_base_path: Path,
                  db_file_name: str = 'ms2_results.sql',
                  query: str = None,
                  trim_spectrum_id: bool = False):

        db_paths = [*project_base_path.rglob(db_file_name)]

        df_list = []

        for db_path in db_paths:

            if query:

                df_tmp = DataLoader.sql_to_df(db_file=db_path, query=query)

            else:

                df_tmp = DataLoader.sql_to_df(db_file=db_path)

            if trim_spectrum_id:

                trim_spectra = lambda x: int(x.split('=')[-1])

                df_tmp['spectrum_id'] = df_tmp['spectrum_id']\
                                        .apply(trim_spectra)

            df_tmp['sample'] = db_path.parts[-4]

            df_tmp['antigen'] = db_path.parts[-3]

            df_tmp['antibody'] = db_path.parts[-2]

            df_list.append(df_tmp)

        self.df = pd.concat(df_list, ignore_index=True)

        del df_list

    def get_dataframe(self) -> DataFrame:

        return self.df

    @classmethod
    def sql_to_df(self,
                  db_file: Path,
                  query: str = 'SELECT * FROM MS2data') -> DataFrame:

        con = sqlite3.connect(db_file)

        return pd.read_sql(query, con)


def load_scores(output_path: Path) -> pd.DataFrame:

    def get_xls(model: str, score_file: Path) -> List:

        xl_file = score_file.parent / f'{model}_xls.txt'

        with open(xl_file, 'r') as f:

            reader = csv.reader(f)

            xl_list = [xl for row in reader for xl in row]

        return xl_list

    score_files = output_path.rglob('*scores.csv')

    df_list = []

    for score_file in score_files:

        df_tmp = pd.read_csv(score_file)

        df_tmp['antibody'] = score_file.parts[-3]

        df_tmp['antigen'] = score_file.parts[-4]

        df_tmp['sample'] = score_file.parts[-5]

        df_tmp['xls'] = df_tmp['model'].apply(get_xls, args=[score_file])

        df_list.append(df_tmp)

    return pd.concat(df_list, ignore_index=True)

def split_peptide(df, antigen_dict):
    # Create boolean mask for NaN values in the 'xls' column
    nan_mask = df['xls'].isna()

    nan = float('nan')

    # Fill NaN values with placeholders for non-peptide rows
    df['ab_pep'] = df['xls'].where(~nan_mask, nan)

    df['ab_pep_xl_index'] = df['xls'].where(~nan_mask, nan)

    df['ag_pep'] = df['xls'].where(~nan_mask, nan)

    df['ag_pep_xl_index'] = df['xls'].where(~nan_mask, nan)

    df['ag_pep_xl_pos'] = df['xls'].where(~nan_mask, nan)

    df['ag_pep_xl_pos_real'] = df['xls'].where(~nan_mask, nan)

    # Only process non-NaN peptides
    non_nan_peptides = df.loc[~nan_mask, 'xls']

    # Split peptides on '--'
    peptides_split = non_nan_peptides.str.split('--', expand=True)

    # Clean ab_pep and ag_pep
    df.loc[~nan_mask, 'ab_pep'] = peptides_split[0].str.strip('-.(0123456789)')

    df.loc[~nan_mask, 'ag_pep'] = peptides_split[1].str.strip('-.(0123456789)')

    # Extract crosslink indices and convert to integers
    ab_pep_xl_index = peptides_split[0].str\
                                    .extract(r'\((\d+)\)')[0]\
                                    .astype(int) - 1

    ag_pep_xl_index = peptides_split[1].str\
                                    .extract(r'\((\d+)\)')[0]\
                                    .astype(int) - 1

    df.loc[~nan_mask, 'ab_pep_xl_index'] = ab_pep_xl_index

    df.loc[~nan_mask, 'ag_pep_xl_index'] = ag_pep_xl_index

    # Compute ag_pep_xl_pos by finding the position of ag_pep in the
    # antigen sequence
    search_antigen = lambda x: [match.start()
                                for match in re.finditer(x[0], antigen_dict[x[1]])]

    df.loc[~nan_mask,
        'ag_pep_xl_pos'] = df.loc[~nan_mask,
                                    ['ag_pep', 'antigen']].apply(lambda x: search_antigen(x), axis=1)

    df = df.explode('ag_pep_xl_pos').reset_index(drop=True)

    df['ag_pep_xl_pos'] += df['ag_pep_xl_index']

    df['ag_pep_xl_pos_real'] = df['ag_pep_xl_pos'] + 1

    return df


@click.command()
@click.version_option(version='1.0a')
@click.option('--pdb_a',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='Antibody (target) PDB')
@click.option('--pdb_b',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='Antigen (binder) PDB')
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
def run_scoring(pdb_a: Path,
              pdb_b: Path,
              output_dir: Path,
              top_xls_file: Path,
              dock_file: Path,
              n_models: int = 10,
              n_filters: int = 3, cut_off: int = 32):
    """
    """
    pm.finish_launching()

    parser = PDBParser(QUIET=True)

    df_scores = load_scores(output_dir)

    lambda_func = lambda x: str(int(x.split('_')[-1])+1)

    df_scores['antibody_model'] = df_scores['antibody'] + '_' +\
                                  df_scores['model'].transform(lambda_func)

    df_scores = df_scores.loc[df_scores['xl_hits'] > 0]

    df_scores_exp = df_scores.explode('xls').reset_index(drop=True).copy()

    antigen_dict = {pdb_b.stem: str(record.seq)
                    for record in SeqIO.parse(str(pdb_b), 'pdb-atom')}

    df_scores_exp = split_peptide(df_scores_exp, antigen_dict)


if __name__ == "__main__":

    run_scoring()
