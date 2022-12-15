#!/usr/bin/env python3
"""Cheetah-MS MS2 analysis.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Remove hardcoded default parameters (delta, intensity, ptm, xlinker)


from pathlib import Path
from sqlalchemy import create_engine

import click
import pandas

from ms2_utils import dda_clean


# xl_file should potentially be changed to a top XL file...
def taxlink(xl_file: Path,
            mgf_file: Path,
            delta: float,
            intensity_filter: float, xlinker_type: int, ptm_type: str) -> Path:
    """...

    ...

    Parameters
    ----------

    """

    # Disuccinimidyl suberate:
    DSS_mass = 138.06808  # mass.calculate_mass(formula='C16H20N2O8')

    # Disuccinimidyl glutarate:
    DSG_mass = 96.02059

    # Ethylene glycol bis(succinimidyl succinate):
    EGS_mass = 226.04774

    xlinker_dic = {1: DSS_mass,
                   2: DSG_mass,
                   3: EGS_mass}

    xlinker_mass = xlinker_dic[xlinker_type]
    fig_threshold = 10
    final_score = 0
    precursor_delta = delta

    cleaned_mgf_file = DDA_clean(xl_file,
                                 mgf_file,
                                 precursor_delta, xlinker_mass, ptm_type)


@click.command()
@click.version_option(version='1.0a')
@click.option('--mgf_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--xl_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--output_dir',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
def run_ms2_analysis(mgf_file: Path, xl_file: Path, output_dir: Path) -> Path:
    """...

    ...

    Parameters
    ----------

    """
    delta = 0.01  # 0.01 or 0.05.

    intensity = 0.0

    x_linker = 1

    ptm = "1"

    taxlink(xl_file=xl_file,
            mgf_file=mgf_file,
            delta=delta,
            intensity_filter=intensity, xlinker_type=x_linker, ptm_type=ptm)

    engine = create_engine('sqlite:///ms2_results.sql')

    query = 'SELECT DISTINCT XL FROM MS2Data WHERE count>10 ORDER BY count DESC'

    ms2_xls_df = pandas.read_sql_query(query, engine)

    with open('top_xls.txt', 'a') as f:

        f.write(ms2_xls_df.to_csv(sep=' ', index=False, header=False))


if __name__ == '__main__':

    run_ms2_analysis()
