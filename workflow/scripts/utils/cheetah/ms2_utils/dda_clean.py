from pathlib import Path
from typing import List

from pyteomics import mgf

from .mgf_utils import fragment_generator


def mgf_writer(mgf_output_file: Path, spectra: dict):
    """Write MGF file.

    ...

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.

    ...

    Parameters
    ----------

    """

    f = mgf_output_file

    title = str(spectra['params']['title'])

    pepmass = str(spectra['params']['pepmass'][0])

    rtinseconds = str(spectra['params']['rtinseconds'])

    charge = (str(spectra["params"]["charge"][0])).replace("+", "")

    print('BEGIN IONS',
            f'TITLE={title}',
            f'PEPMASS={pepmass}',
            f'RTINSECONDS={rtinseconds}',
            f'CHARGE={charge}', sep='\n', file=f)

    for i in range(len(spectra['m/z array'])):

        mz = str(spectra["m/z array"][i])

        intensity = str(spectra["intensity array"][i])

        print(f'{mz}', f'{intensity}', sep=' ', file=f)

    print('END IONS\n', file=f)

def DDA_clean(xl_list: List[str],
              mgf_file: Path,
              precursor_delta: float, xlinker_mass: int, ptm_type: str) -> Path:
    """Clean DDA MGF file.

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.

    ...

    Parameters
    ----------
    xl_list : list of strings
        List strings with Kojak formatted cross-links.
    mgf_file : pathlib.Path object
        Path to MGF with MS2 data.
    precursor_delta : float
    xlinker_mass : int
    ptm_type : str
        Valid input {'1'..'12'}.

    Returns
    -------
    pathlib.Path object
        Path to cleaned MGF file.
    """

    xls = xl_list

    # Read the MGF (MS/MS) file.
    mgf_dict_list = list(mgf.read(str(mgf_file)))

    output_file = mgf_file.parent / f'{mgf_file.stem}_cleaned.mgf'

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
