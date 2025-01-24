#!/usr/bin/env python3

from pathlib import Path
from typing import List

import pandas as pd
import pymzml
from pandas import Series

from . import fragment_generator


def mgf_writer(mgf_output_file: Path,
               series: Series) -> None:
    """Writes a single MS/MS spectrum to an MGF file.

    This function takes a dictionary containing the metadata and peak
    information of a single MS/MS spectrum and writes it to the specified MGF
    file.

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.

    Args:
        mgf_output_file (Path): Path to the output MGF file (opened in append mode).
        spectra (dict): Dictionary containing the MS/MS spectrum data.
            - 'params': dictionary with spectrum metadata (title, pepmass, rtinseconds, charge).
            - 'm/z array': list of m/z values for the spectrum peaks.
            - 'intensity array': list of intensity values for the spectrum peaks.
    """
    title = series['spectra_id']

    pepmass = series['mz']

    pepintensity = series['i']

    rtinseconds = series['rt']

    charge = series['charge']

    mgf_output_file.write('BEGIN IONS\n')

    mgf_output_file.write(f'TITLE={title}\n')

    mgf_output_file.write(f'PEPMASS={pepmass}\n')

    mgf_output_file.write(f'PEPINTENSITY={pepintensity}\n')

    mgf_output_file.write(f'RTINSECONDS={rtinseconds}\n')

    mgf_output_file.write(f'CHARGE={charge}\n')

    for mz, intensity in series['ms2']:

        mgf_output_file.write(f'{mz} {intensity}\n')

    mgf_output_file.write('END IONS\n')

def dda_filter(xl_list: List[str],
               mzml_file: Path,
               output_dir: Path,
               precursor_delta: float,
               xlinker_mass: int, xlinker: int, ptm_type: str) -> Path:
    """Filters an MGF file based on precursor masses matching cross-links.

    This function takes a list of potential cross-links (XLs), an MGF file containing MS/MS spectra,
    an output directory, a mass tolerance (delta) for precursor ion matching, the XL linker mass,
    and the PTM type (currently supports unmodified peptides only). It performs the following steps:

    1. Reads the MGF file to access individual spectra.
    2. Iterates through each XL in the list.
        - Generates theoretical fragment ions for the XL sequence.
    3. Iterates through each spectrum in the MGF file.
        - Compares the precursor mass of the spectrum to the theoretical precursor masses
          of the light and heavy forms of the XL (considering the linker mass) within the specified tolerance.
        - If a match is found, the entire spectrum is written to a new filtered MGF file.

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.

    Args:
        xl_list (List[str]): List of strings containing Kojak-formatted cross-links.
        mgf_file (Path): Path to the MGF file containing MS/MS spectra.
        output_dir (Path): Path to the output directory for storing the filtered MGF file.
        precursor_delta (float): Mass tolerance (delta) for precursor ion matching.
        xlinker_mass (int): Mass of the XL linker used.
        ptm_type (str): PTM type considered (currently supports unmodified peptides only, "1").

    Returns:
        Path: Path to the generated filtered MGF file.
    """
    xls = xl_list

    # Read the MGF (MS/MS) file.
    mzml = pymzml.run.Reader(str(mzml_file))

    id_list = []

    rt_list = []

    ms2_list = []

    mz_list = []

    i_list = []

    charge_list = []

    for spectra in mzml:

        if spectra.ms_level == 2:

            for precursor in spectra.selected_precursors:

                id_list.append(spectra.ID)

                ms2_list.append(spectra.peaks('raw'))

                rt_list.append(spectra.scan_time_in_minutes() * 60)

                mz_list.append(precursor['mz'])

                i_list.append(precursor['i'])

                charge_list.append(precursor['charge'])

    df = pd.DataFrame({'spectra_id': id_list,
                       'mz': mz_list,
                       'i': i_list,
                       'charge': charge_list,
                       'rt': rt_list, 'ms2': ms2_list})

    output_file = output_dir / f'{mzml_file.stem}_filtered.mzML'

    index_memory = []

    with open(output_file, 'w') as f:

        for num_xl, xl in enumerate(xls):

            (precursor_dict,
             fragment_all,
             mz_light_all,
             mz_heavy_all,
             p1, p2) = fragment_generator.fragment_generator(xl,
                                          xlinker_mass,
                                          xlinker,
                                          ptm_type)

            for charge, mass_list in precursor_dict.items():

                index = []

                for mass in mass_list:

                    index.extend(df.loc[(df['charge'] == charge) &
                                        (df['mz'].between(mass-precursor_delta,
                                                          mass+precursor_delta,
                                                          'neither'))]\
                                .index)

                for df_i in index:

                    if df_i not in index_memory:

                        mgf_writer(mgf_output_file=f, series=df.loc[df_i])

                        index_memory.append(df_i)

    return output_file
