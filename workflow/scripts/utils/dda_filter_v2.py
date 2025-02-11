#!/usr/bin/env python3

import numpy as np
from pathlib import Path
from typing import List

import pandas as pd
import pymzml
from pandas import Series

from . import fragment_generator


def check_floats_in_integer_ranges(float_series, range_series, offset=0.01):
    """
    Checks if a Series of floats falls within ranges defined by integers in
    another Series (integer ± offset).

    Args:
        float_series: Pandas Series of floats to be checked.
        range_series: Pandas Series of integers defining the center of ranges.
        offset: The offset around each integer to define the range.

    Returns:
        Pandas Series of booleans.
    """

    def is_in_range(value):
        for center in range_series:
            if center - offset <= value <= center + offset:
                return True
        return False

    return float_series[float_series.apply(is_in_range)].index

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

    ms2_array = np.array(series['ms2'])

    np.savetxt(mgf_output_file, ms2_array, delimiter=' ', fmt='%15.30f')

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

    data = []

    for spectra in mzml:

        if spectra.ms_level == 2:

            for precursor in spectra.selected_precursors:

                data.append({'spectra_id': spectra.ID,
                             'mz': precursor['mz'],
                             'i': precursor['i'],
                             'charge': precursor['charge'],
                             'rt': spectra.scan_time_in_minutes() * 60,
                             'ms2': np.array(spectra.peaks('raw'))})

    df = pd.DataFrame(data)

    output_file = output_dir / f'{mzml_file.stem}_filtered.mgf'

    xls_precursor_dict = {i: [] for i in range(3, 9)}

    for xl in xls:

        (precursor_dict,
            *_) = fragment_generator.fragment_generator(xl,
                                                        xlinker_mass,
                                                        xlinker, ptm_type)

        for charge, mass_list in precursor_dict.items():

            xls_precursor_dict[charge].extend(mass_list)

    filter_index = []

    for charge in range(3, 9):

        float_series = df.loc[df['charge'] == charge]['mz']

        range_series = xls_precursor_dict[charge]

        index = check_floats_in_integer_ranges(float_series,
                                               range_series, precursor_delta)

        filter_index.extend(index.to_list())

    with open(output_file, 'w') as f:

        for i in filter_index:

            mgf_writer(f, df.loc[i])

    return output_file
