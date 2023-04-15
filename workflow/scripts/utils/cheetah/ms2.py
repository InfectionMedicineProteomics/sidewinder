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

from ms2_utils import dda_filter
from ms2_utils import fig_maker
from ms2_utils.mgf_utils import fragment_generator


# xl_file should potentially be changed to a top XL file...
def taxlink(all_xls_file: Path,
            mgf_file: Path,
            output_dir: Path,
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

    xlinker_dict = {1: DSS_mass,
                   2: DSG_mass,
                   3: EGS_mass}

    xlinker_mass = xlinker_dict[xlinker_type]

    fig_threshold = 10  # Set number of figures.

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

    con = sqlite3.connect(sql_db_file)

    c = con.cursor()

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

                    # mz_difference_Light = 1.000

                    # mz_difference_Heavy = 1.000

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

                                    # Intensity-based scoring of the detected peaks.
                                    if spectra['intensity array'][y] >= maxintensity*3/4:

                                        mgf_spec_score[spec_idx] += 12

                                    elif spectra['intensity array'][y] >= maxintensity/2:

                                        mgf_spec_score[spec_idx] += 8

                                    elif spectra['intensity array'][y] >= maxintensity/4:

                                        mgf_spec_score[spec_idx] += 4

                                    elif spectra['intensity array'][y]>=intensity_filter:

                                        mgf_spec_score[spec_idx] += 2

                    elif (mz_difference_Heavy <= delta):

                        print(spectra['params']['title'], file=f)

                        for num_mz, anymZ in enumerate(mz_heavy_all):

                            for y in range(len(spectra['m/z array'])):
                                min_MZ = abs(anymZ - spectra['m/z array'][y])
                                if min_MZ <= delta:
                                    maxintensity = max(spectra['intensity array'][0:])
                                    mgf_spec_list[num_mz] = spectra

                                    # Intensity-based scoring of the detected peaks.
                                    if spectra['intensity array'][y] >= maxintensity*3/4:
                                        mgf_spec_score[spec_idx] += 12
                                    elif spectra['intensity array'][y] >= maxintensity/2:
                                        mgf_spec_score[spec_idx] += 8
                                    elif spectra['intensity array'][y] >= maxintensity/4:
                                        mgf_spec_score[spec_idx] += 4
                                    elif spectra['intensity array'][y]>=intensity_filter:
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

                        elif (mz_difference_Heavy <= delta): # we fixed delta for precursor to 0.01

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
                            ,(xl, str(all_xls_file), spectrum_id, spectrum_num, delta, spec_charge, heavy_light_flag, \
                            max(mgf_spec_score), coverage, ",".join(str(x) for x in covered_Frags), \
                            ",".join(str(w) for w in covered_mz), ",".join(str(x) for x in covered_intensity), \
                            ",".join(str(x) for x in main_spectra), ",".join(str(x) for x in main_intensity), \
                            len(covered_Frags) ) )

                            con.commit()

                        if len(covered_mz) >= fig_threshold:

                            fig_maker.fig_maker(main_spectra,
                                                main_intensity,
                                                covered_Frags,
                                                covered_mz,
                                                p1,
                                                p2,
                                                xl, num_mgf, delta, output_dir)

    con.close()

    return sql_db_file

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

    sql_db_path = taxlink(all_xls_file=xl_file,
                          mgf_file=mgf_file,
                          output_dir=output_dir,
                          delta=delta,
                          intensity_filter=intensity,
                          xlinker_type=x_linker, ptm_type=ptm)

    engine = create_engine(f'sqlite:///{str(sql_db_path)}')

    # Save the XLs that have the most support from the MS2 spectra.
    query = 'SELECT DISTINCT XL FROM MS2Data WHERE count>10 ORDER BY count DESC'

    ms2_xls_df = pandas.read_sql_query(query, engine)

    top_xls_file = output_dir / 'top_xls.txt'

    with open(top_xls_file, 'w') as f:

        f.write(ms2_xls_df.to_csv(sep=' ', index=False, header=False))


if __name__ == '__main__':

    run_ms2_analysis()
