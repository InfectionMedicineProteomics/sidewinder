#!/usr/bin/env python3
"""Cheetah-MS MS2 analysis.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Remove hardcoded default parameters (delta, intensity, ptm, xlinker)
# - Clean up the taxlink function


import sqlite3
from pathlib import Path
from sqlalchemy import create_engine

import click
import pandas
from pyteomics import mgf

from utils import dda_filter
from utils import fig_maker
from utils import fragment_generator


# xl_file should potentially be changed to a top XL file...
def taxlink(all_xls_file: Path,
            mgf_file: Path,
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
    # Disuccinimidyl suberate:
    DSS_mass = 138.06808  # mass.calculate_mass(formula='C16H20N2O8')

    # Disuccinimidyl glutarate:
    DSG_mass = 96.02059

    # Ethylene glycol bis(succinimidyl succinate):
    EGS_mass = 226.04774

    xlinker_dict = {1: DSS_mass,
                   2: DSG_mass,
                   3: EGS_mass}

    xlinker_mass = xlinker_dict[xlinker_type]

    fig_threshold = 10  # Min n spectra to generate fig.

    # Importing XLs file.
    with open(all_xls_file, 'r') as f:

        all_xls_list = f.read().splitlines()  # Row -> element of the list.

    # Considering the list of input XLs, we remove all extra spectra
    # from input MGF according to the pre-cursor m/z and make a new
    # filtered version to analyze.
    filtered_mgf_file = dda_filter.DDA_filter(all_xls_list,
                                             mgf_file,
                                             output_dir,
                                             delta,
                                             xlinker_mass, ptm_type)

    # Build SQLite3 database.
    sql_db_file = output_dir / 'ms2_results.sql'

    conn = sqlite3.connect(sql_db_file)

    c = conn.cursor()

    # Create table:
    c.execute('''CREATE TABLE IF NOT EXISTS MS2Data
                 (XL text, mgf_file text, spectrum_id text, spectrum_num integer, \
                 delta real, pre_charge integer, H_L text, fragSc integer, \
                 coverage real, covered_Frags text, covered_Mz text, covered_int text,\
                 main_Mz text, main_int text, count integer)''')

    c.execute('PRAGMA journal_mode = WAL')

    c.execute('PRAGMA synchronous = NORMAL')

    # Read in cleaned MGF file.
    mgf_dict_list = list(mgf.read(str(filtered_mgf_file)))

    output_file = output_dir / 'detected_spectra.txt'

    spec_dir = output_dir / 'top_spectra'

    spec_dir.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:

        for num_xl, xl in enumerate(all_xls_list):

            (precursor_dict,
             fragment_all,
             mz_light_all,
             mz_heavy_all,
             p1,
             p2) = fragment_generator.fragment_generator(xl,
                                                         xlinker_mass, ptm_type)

            mgf_spec_score = [0] * len(mgf_dict_list)

            ## Considering all fragments and find the match peaks in MS/MS data
            mgf_spec_list = [None] * len(mz_light_all)

            print(xl, file=f)

            for spec_idx, spectra in enumerate(mgf_dict_list):

                spec_pepmass = spectra['params']['pepmass'][0]

                spec_charge = spectra['params']['charge'][0]

                if spec_charge in range (3,9):

                    mz_difference_Light = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][0:4])

                    mz_difference_Heavy = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][5:9])

                    if (mz_difference_Light <= delta):

                        print(spectra['params']['title'], file=f)

                        for num_mz, anymZ in enumerate(mz_light_all):

                            for y in range(len(spectra['m/z array'])):

                                min_MZ = abs(anymZ - spectra['m/z array'][y])

                                if min_MZ <= delta:

                                    maxintensity = max(spectra['intensity array'][0:])

                                    mgf_spec_list[num_mz] = spectra

                                    # Intensity-based scoring
                                    # of the detected peaks:
                                    intensity = spectra['intensity array'][y]

                                    if intensity >= maxintensity*3/4:

                                        mgf_spec_score[spec_idx] += 12

                                    elif intensity >= maxintensity/2:

                                        mgf_spec_score[spec_idx] += 8

                                    elif intensity >= maxintensity/4:

                                        mgf_spec_score[spec_idx] += 4

                                    elif intensity >= intensity_filter:

                                        mgf_spec_score[spec_idx] += 2

                    elif (mz_difference_Heavy <= delta):

                        print(spectra['params']['title'], file=f)

                        for num_mz, anymZ in enumerate(mz_heavy_all):

                            for y in range(len(spectra['m/z array'])):

                                min_MZ = abs(anymZ - spectra['m/z array'][y])

                                if min_MZ <= delta:

                                    maxintensity = max(spectra['intensity array'][0:])

                                    mgf_spec_list[num_mz] = spectra

                                    # Intensity-based scoring
                                    # of the detected peaks:
                                    intensity = spectra['intensity array'][y]

                                    if intensity >= maxintensity*3/4:

                                        mgf_spec_score[spec_idx] += 12

                                    elif intensity >= maxintensity/2:

                                        mgf_spec_score[spec_idx] += 8

                                    elif intensity >= maxintensity/4:

                                        mgf_spec_score[spec_idx] += 4

                                    elif intensity >= intensity_filter:

                                        mgf_spec_score[spec_idx] += 2

            spectrum_id = "NA"

            spectrum_num = 0

            peaks_num = 0

            coverage = 0.0

            for num_mgf, item_mgf in enumerate(mgf_dict_list):

                if max(mgf_spec_score) != 0:

                    if mgf_spec_score[num_mgf] == max(mgf_spec_score):

                        heavy_light_flag = ''

                        spec_pepmass = item_mgf['params']['pepmass'][0]

                        spec_charge = item_mgf['params']['charge'][0]

                        # mz_difference_Light = 1.000

                        # mz_difference_Heavy = 1.000

                        mz_difference_Light = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][0:4])

                        mz_difference_Heavy = min(abs(spec_pepmass - pre_mz) for pre_mz in precursor_dict[spec_charge][5:9])

                        spectrum_id = item_mgf['params']['title']

                        spectrum_num = num_mgf

                        peaks_num = len(item_mgf['m/z array'])

                        ## To send to fig_maker:
                        main_spectra = item_mgf['m/z array']

                        main_intensity = item_mgf['intensity array']

                        covered_Frags = []

                        covered_mz = []

                        covered_intensity = []

                        if (mz_difference_Light <= delta): # we fixed delta for precursor to 0.01

                            heavy_light_flag = 'Light'

                            for numz, mz in enumerate(mz_light_all):

                                min_MZ2 = min(abs(mz - y) for y in item_mgf['m/z array'][0:])

                                if min_MZ2 <= delta:

                                    covered_Frags.append(fragment_all[numz])

                                    covered_mz.append(mz_light_all[numz])

                        elif (mz_difference_Heavy <= delta):

                            heavy_light_flag = 'Heavy'

                            for numz, mz in enumerate(mz_heavy_all):

                                min_MZ2 = min(abs(mz - y) for y in item_mgf['m/z array'][0:])

                                if min_MZ2 <= delta:

                                    covered_Frags.append(fragment_all[numz])

                                    covered_mz.append(mz_heavy_all[numz])

                        for y_num, ymz in enumerate(covered_mz):

                            for mz_num,anymz in enumerate(main_spectra):

                                min_MZ = abs(anymz - ymz)

                                if min_MZ <= delta:

                                    if main_intensity[mz_num] != 0.0:

                                        covered_intensity.append(main_intensity[mz_num])
                                        break

                        if peaks_num != 0:

                            coverage = (len(covered_Frags) * 100)/peaks_num

                        if (spectrum_id != 'NA' and len(covered_mz) >= 5):

                            c.execute("INSERT INTO MS2Data(XL, mgf_file,spectrum_id, spectrum_num, \
                            delta, pre_charge, H_L, fragSc, coverage, covered_Frags, covered_Mz, covered_int, \
                            main_Mz, main_int, count) VALUES \
                            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)" \
                            ,(xl, str(filtered_mgf_file), spectrum_id, spectrum_num, delta, spec_charge, heavy_light_flag, \
                            max(mgf_spec_score), coverage, ",".join(str(x) for x in covered_Frags), \
                            ",".join(str(w) for w in covered_mz), ",".join(str(x) for x in covered_intensity), \
                            ",".join(str(x) for x in main_spectra), ",".join(str(x) for x in main_intensity), \
                            len(covered_Frags) ) )

                            conn.commit()

                        if len(covered_mz) >= fig_threshold:

                            fig_maker.fig_maker(main_spectra,
                                                main_intensity,
                                                covered_Frags,
                                                covered_mz,
                                                p1,
                                                p2,
                                                xl, num_mgf, delta, spec_dir)

    return conn

@click.command()
@click.version_option(version='1.0a')
@click.option('--mgf_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--xl_file',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--x_linker',
              required=False, default=1, type=int,
              help='Options: 1=DSS, 2=DSG, 3=EGS')
@click.option('--mass_delta_cutoff', type=float, default=0.01, help='')
@click.option('--output_dir',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
def run_ms2_analysis(mgf_file: Path,
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
        mgf_file (Path): Path to the MGF file containing MS/MS spectra.
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
                          mgf_file=mgf_file,
                          output_dir=output_dir,
                          delta=delta,
                          intensity_filter=intensity,
                          xlinker_type=x_linker, ptm_type=ptm)

    # Save the XLs that have the most support from the MS2 spectra.
    query = ('SELECT DISTINCT XL FROM MS2Data '
             'WHERE count>=10 ORDER BY count DESC')

    ms2_xls_df = pandas.read_sql_query(query, sql_db_conn)

    sql_db_conn.close()

    top_xls_file = output_dir / 'top_xls.txt'

    with open(top_xls_file, 'w') as f:

        f.write(ms2_xls_df.to_csv(sep=' ', index=False, header=False))


if __name__ == '__main__':

    run_ms2_analysis()
