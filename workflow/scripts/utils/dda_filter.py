#!/usr/bin/env python3

from pathlib import Path
from typing import List

from pyteomics import mgf

from . import fragment_generator


def mgf_writer(mgf_output_file: Path, spectra: dict) -> None:
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
    with open(mgf_output_file, 'a') as f:  # Append mode for multiple spectra.

        f.write('BEGIN IONS\n')

        f.write(f'TITLE={spectra["params"]["title"]}\n')

        f.write(f'PEPMASS={spectra["params"]["pepmass"][0]}\n')

        f.write(f'RTINSECONDS={spectra["params"]["rtinseconds"]}\n')

        f.write(f'CHARGE={(str(spectra["params"]["charge"][0]))\
                          .replace("+", "")}\n')

        for i in range(len(spectra['m/z array'])):

            mz = str(spectra["m/z array"][i])

            intensity = str(spectra["intensity array"][i])

            f.write(f'{mz} {intensity}\n')

        f.write('END IONS\n\n')

def DDA_filter(xl_list: List[str],
              mgf_file: Path,
              output_dir: Path,
              precursor_delta: float, xlinker_mass: int, ptm_type: str) -> Path:
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
    mgf_dict_list = list(mgf.read(str(mgf_file)))

    output_file = output_dir / f'{mgf_file.stem}_filtered.mgf'

    with open(output_file, 'w') as f:

        for num_xl, xl in enumerate(xls):

            (precursor_dict,
             fragment_all,
             mz_light_all,
             mz_heavy_all,
             p1, p2) = fragment_generator.fragment_generator(xl,
                                                             xlinker_mass,
                                                             ptm_type)

            for spectra in mgf_dict_list:  # Removed enumeration.

                try:

                    spec_pepmass = spectra['params']['pepmass'][0]

                    spec_charge = spectra['params']['charge'][0]

                    if spec_charge in range(3, 9):

                        # I don't understand why these are declared...
                        #mz_diff_light = 1.000

                        #mz_diff_heavy = 1.000

                        mz_diff_light = min(abs(spec_pepmass - pre_mz)
                                            for pre_mz
                                            in precursor_dict[spec_charge][0:4])

                        mz_diff_heavy = min(abs(spec_pepmass - pre_mz)
                                            for pre_mz
                                            in precursor_dict[spec_charge][5:9])

                        if mz_diff_light <= precursor_delta:

                            mgf_writer(mgf_output_file=f, spectra=spectra)

                        elif mz_diff_heavy <= precursor_delta:

                            mgf_writer(mgf_output_file=f, spectra=spectra)

                except TypeError as error:

                    if "'NoneType' object is not subscriptable" in str(error):

                        print('STRANGE!',
                              'Make sure this is not breaking anything')

                    else:

                        raise error

    return output_file
