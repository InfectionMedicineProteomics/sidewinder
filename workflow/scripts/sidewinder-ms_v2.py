#!/usr/bin/env python3
"""Cheetah-MS MS2 analysis.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Remove hardcoded default parameters (delta, intensity, ptm, xlinker)
# - Clean up the taxlink function


import math
import numpy as np
import sqlite3
from pathlib import Path
from sqlalchemy import create_engine

import click
import pandas as pd
from pyteomics import mgf
from scipy.stats import zscore

from utils import dda_filter_v2
from utils import fig_maker
from utils import fragment_generator


# Constants:
# Disuccinimidyl suberate:
DSS_mass = 138.06808  # mass.calculate_mass(formula='C16H20N2O8')

# Disuccinimidyl glutarate:
DSG_mass = 96.02059

# Ethylene glycol bis(succinimidyl succinate):
EGS_mass = 226.04774

xlinker_dict = {1: DSS_mass, 2: DSG_mass, 3: EGS_mass}

xlinker_name_dict = {1: 'DSS', 2: 'DSG', 3: 'EGS'}


# xl_file should potentially be changed to a top XL file...
def taxlink(all_xls_file: Path,
            mzml_file: Path,
            output_dir: Path,
            delta: float,
            intensity_filter: float, xlinker_type: int, ptm_type: str) -> Path:
    """
    Analyzes MS/MS data (MGF file) to identify spectra supporting potential XL
    interactions (XLs) defined in a separate file.

    This function takes an MGF file containing MS/MS spectra, a list of potential XL links
    (one per line in a text file), an output directory, mass tolerance (delta) for precursor
    and fragment ion matching, an intensity filter for considering matched peaks, an XL linker
    type (DSS, DSG, or EGS), and a PTM type (currently supports unmodified peptides only).

    It performs the following steps:

    1. Imports the list of XLs from the provided file.
    2. Filters the MGF file based on precursor masses matching the XLs within the specified tolerance.
    3. Creates a SQLite database for storing analysis results.
    4. Iterates through each XL:
        - Generates theoretical fragment ions for the XL sequence considering the linker mass and PTM type.
        - Analyzes the filtered MGF data to identify spectra matching both the precursor mass
          and a minimum number of fragment ions for the XL, considering the mass tolerance.
        - Scores the matching spectra based on the number and intensity of matched fragment ions.
        - Stores the results (XL, matched spectrum ID, fragment coverage, etc.) in the SQLite database.
        - Optionally creates figures for top-scoring spectra (with a minimum fragment coverage threshold).
    5. Saves the XLs with the most supported spectra (highest scoring) to a separate text file.

    Args:
        all_xls_file (Path): Path to the text file containing potential XL interactions (one per line).
        mgf_file (Path): Path to the MGF file containing MS/MS spectra.
        output_dir (Path): Path to the output directory for storing results.
        delta (float): Mass tolerance (delta) for precursor and fragment ion matching.
        intensity_filter (float): Minimum intensity filter for considering matched peaks in spectra.
        xlinker_type (int): Type of XL linker used (1: DSS, 2: DSG, 3: EGS).
        ptm_type (str): Type of PTM considered (currently supports unmodified peptides only, "1").

    Returns:
        Path: A handle to the opened SQLite database connection (closed upon script completion).
    """
    xlinker_name = xlinker_name_dict[xlinker_type]

    xlinker_mass = xlinker_dict[xlinker_type]

    fig_threshold = 10  # Min n spectra to generate fig.

    # Importing XLs file.
    with open(all_xls_file, 'r') as f:

        all_xls_list = f.read().splitlines()  # Row -> element of the list.

    # Considering the list of input XLs, we remove all extra spectra
    # from input MGF according to the pre-cursor m/z and make a new
    # filtered version to analyze.
    filtered_mgf_file = dda_filter_v2.dda_filter(all_xls_list,
                                                 mzml_file,
                                                 output_dir,
                                                 delta,
                                                 xlinker_mass,
                                                 xlinker_type, ptm_type)

    # Build SQLite3 database.
    sql_db_file = output_dir / 'ms2_results.sql'

    conn = sqlite3.connect(str(sql_db_file))

    c = conn.cursor()

    # Create table:
    c.execute('''CREATE TABLE IF NOT EXISTS MS2Data
                 (XL text, linker text, XL_mass real, XL_intensity real, mgf_file text, spectrum_id text, spectrum_num integer, \
                 delta real, pre_charge integer, H_L text, fragSc integer, \
                 coverage real, covered_Frags text, covered_Mz text, covered_int text,\
                 main_Mz text, main_int text, count integer, S4b_score real, log_intensity real, S4b_zscore real, log_intensity_zscore real)''')

    c.execute('PRAGMA journal_mode = WAL')

    c.execute('PRAGMA synchronous = NORMAL')

    # Read in cleaned MGF file.
    mgf_dict_list = list(mgf.read(str(filtered_mgf_file)))

    num_spectra = len(mgf_dict_list)

    output_file = output_dir / 'detected_spectra.txt'

    spec_dir = output_dir / 'top_spectra'

    spec_dir.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:

        for xl in all_xls_list:

            (precursor_dict,
             fragment_all,
             mz_light_all,
             mz_heavy_all,
             p1,
             p2) = fragment_generator.fragment_generator(xl,
                                                         xlinker_mass,
                                                         xlinker_type, ptm_type)

            spec_score = np.zeros(num_spectra, dtype=int)

            print(xl, file=f)

            for spec_idx, spectra in enumerate(mgf_dict_list):

                spec_pepmass = spectra['params']['pepmass'][0]

                spec_charge = spectra['params']['charge'][0]

                if 3 <= spec_charge <= 8:

                    mz_difference_Light = (min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][0:5]))

                    mz_difference_Heavy = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][5:9])

                    if (mz_difference_Light <= delta
                        or mz_difference_Heavy <= delta):

                        print(spectra['params']['title'], file=f)

                        mz_to_match = (mz_light_all
                                       if mz_difference_Light <= delta
                                       else mz_heavy_all)

                        for num_mz, anymZ in enumerate(mz_to_match):

                            for y in range(len(spectra['m/z array'])):

                                min_MZ = abs(anymZ - spectra['m/z array'][y])

                                if min_MZ <= delta:

                                    maxintensity = max(spectra['intensity array'])

                                    # Intensity-based scoring
                                    # of the detected peaks:
                                    intensity = spectra['intensity array'][y]

                                    if intensity >= maxintensity*3/4:

                                        spec_score[spec_idx] += 12

                                    elif intensity >= maxintensity/2:

                                        spec_score[spec_idx] += 8

                                    elif intensity >= maxintensity/4:

                                        spec_score[spec_idx] += 4

                                    elif intensity >= intensity_filter:

                                        spec_score[spec_idx] += 2

            max_score_idx = np.argmax(spec_score)

            max_score = spec_score[max_score_idx]

            if max_score > 0:

                spectra = mgf_dict_list[max_score_idx]

                spectrum_id = "NA"

                spectrum_num = 0

                peaks_num = 0

                coverage = 0.0

                spec_pepmass = spectra['params']['pepmass'][0]

                spec_pepi = float(spectra['params']['pepintensity'])

                spec_charge = spectra['params']['charge'][0]

                mz_difference_Light = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][0:5])

                mz_difference_Heavy = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][5:9])

                heavy_light_flag = ('Light'
                                    if mz_difference_Light <= delta
                                    else 'Heavy')

                spectrum_id = spectra['params']['title']

                spectrum_num = max_score_idx

                peaks_num = len(spectra['m/z array'])

                ## To send to fig_maker:
                main_spectra = spectra['m/z array']

                main_intensity = spectra['intensity array']

                mz_to_match = (mz_light_all
                               if heavy_light_flag == 'Light' else mz_heavy_all)

                covered_Frags = []

                covered_mz = []

                covered_intensity = []

                for numz, mz in enumerate(mz_to_match):

                    min_MZ2 = min(abs(mz - y) for y in spectra['m/z array'])

                    if min_MZ2 <= delta:

                        covered_Frags.append(fragment_all[numz])

                        covered_mz.append(mz_to_match[numz])

                for y_num, ymz in enumerate(covered_mz):

                    for mz_num,anymz in enumerate(main_spectra):

                        min_MZ = abs(anymz - ymz)

                        if min_MZ <= delta:

                            if main_intensity[mz_num] != 0.0:

                                covered_intensity.append(main_intensity[mz_num])

                                break

                ms2_ions = len(covered_Frags)

                if peaks_num != 0:

                    coverage = (ms2_ions * 100) / peaks_num

                if (spectrum_id != 'NA' and len(covered_mz) >= 5):

                    i_max = max(main_intensity)

                    i_norm_sum = sum([float(i) / i_max
                                        for i in covered_intensity])

                    S4b_score = (ms2_ions+i_norm_sum)/(len(p1)+len(p2))

                    log_intensity = math.log10(spec_pepi + 1)

                    c.execute("INSERT INTO MS2Data(XL, \
                                                    linker, \
                                                    XL_mass, \
                                                    XL_intensity, \
                                                    mgf_file, \
                                                    spectrum_id, \
                                                    spectrum_num, \
                                                    delta, \
                                                    pre_charge, \
                                                    H_L, \
                                                    fragSc, \
                                                    coverage, \
                                                    covered_Frags, \
                                                    covered_Mz, \
                                                    covered_int, \
                                                    main_Mz, \
                                                    main_int, \
                                                    count, \
                                                    S4b_score, \
                                                    log_intensity) \
                                VALUES(?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, \
                                        ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                (xl,
                                xlinker_name,
                                spec_pepmass,
                                spec_pepi,
                                str(filtered_mgf_file),
                                spectrum_id,
                                spectrum_num,
                                delta,
                                spec_charge,
                                heavy_light_flag,
                                max(spec_score),
                                coverage,
                                ",".join(str(x) for x in covered_Frags),
                                ",".join(str(w) for w in covered_mz),
                                ",".join(str(x)
                                        for x in covered_intensity),
                                ",".join(str(x) for x in main_spectra),
                                ",".join(str(x)
                                        for x in main_intensity),
                                len(covered_Frags),
                                S4b_score, log_intensity))

                    conn.commit()

                if len(covered_mz) >= fig_threshold:

                    fig_maker.fig_maker(main_spectra,
                                        main_intensity,
                                        covered_Frags,
                                        covered_mz,
                                        p1,
                                        p2,
                                        xl, max_score_idx, delta, spec_dir)

    return conn

@click.command()
@click.version_option(version='1.0a')
@click.option('--mzml_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--xl_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--x_linker',
              required=False,
              default=1, type=int, help='Options: 1=DSS, 2=DSG, 3=EGS')
@click.option('--mass_delta_cutoff', type=float, default=0.01, help='')
@click.option('--output_dir',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
def run_ms2_analysis(mzml_file: Path,
                     xl_file: Path,
                     x_linker: int,
                     mass_delta_cutoff: float, output_dir: Path) -> Path:
    """
    Runs the Cheetah-MS MS2 analysis.

    This function serves as the main entry point for the script, and parses
    command-line arguments, calls the `taxlink` function to perform the main
    analysis, and saves the identified top-scoring XLs (with the most
    supporting spectra) to a separate file.

    Args:
        mzml_file (Path): Path to the mzML file containing MS/MS spectra.
        xl_file (Path): Path to the text file containing potential XL links.
        x_linker (int): XL linker type (1: DSS, 2: DSG, 3: EGS).
        mass_delta_cutoff (float): Mass tolerance (delta) for precursor and
            fragment ion matching.
        output_dir (Path): Path to the output directory for storing results.

    Returns:
        Path: A handle to the opened SQLite database connection (closed upon script completion).
    """
    delta = mass_delta_cutoff  # 0.01 or 0.05.

    intensity = 0.0

    ptm = "1"

    sql_db_conn = taxlink(all_xls_file=xl_file,
                          mzml_file=mzml_file,
                          output_dir=output_dir,
                          delta=delta,
                          intensity_filter=intensity,
                          xlinker_type=x_linker, ptm_type=ptm)

    # calculate the zscores! and add back to the db
    c = sql_db_conn.cursor()

    selection_query = ('SELECT S4b_score, log_intensity FROM MS2Data')

    s4b_log_i_df = pd.read_sql_query(selection_query, sql_db_conn)

    s4b_zscore = zscore(s4b_log_i_df['S4b_score'])

    log_i_zscore = zscore(s4b_log_i_df['log_intensity'])

    # insert_tuple_list = [(x, y) for x, y in zip(s4b_zscore, log_i_zscore)]

    for i, insert_tuple in enumerate(zip(s4b_zscore, log_i_zscore)):

        rowid = i + 1

        c.execute(f'UPDATE MS2Data SET S4b_zscore=?, log_intensity_zscore=? '
                  f'WHERE ROWID={rowid}', insert_tuple)

    sql_db_conn.commit()

    filter_by = f'S4b_zscore >= 0'  # AND log_intensity_zscore >= 0
                                    # - to also filter on MS1
                                    # - don't think this is good for XL-MS

    # Save the XLs that have the most support from the MS2 spectra.
    filter_query = (f'SELECT DISTINCT XL FROM MS2Data '
                    f'WHERE {filter_by}')

    ms2_xls_df = pd.read_sql_query(filter_query, sql_db_conn)

    sql_db_conn.close()

    top_xls_file = output_dir / 'top_xls.txt'

    with open(top_xls_file, 'w') as f:

        f.write(ms2_xls_df.to_csv(sep=' ', index=False, header=False))


if __name__ == '__main__':

    run_ms2_analysis()
