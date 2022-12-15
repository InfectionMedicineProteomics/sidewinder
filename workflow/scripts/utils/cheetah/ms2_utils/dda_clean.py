from pathlib import Path

from pyteomics import mgf

from mgf_utils import fragment_generator


def mgf_writer(mgf_output_file: Path, spectra: dict):
    """Write MGF file.

    Keyword inputs:
    mgf_output_file -- output file path
    spectra -- MGF spectra dict

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.
    """

    title = str(spectra['params']['title'])

    pepmass = str(spectra['params']['pepmass'][0])

    rtinseconds = str(spectra['params']['rtinseconds'])

    charge = (str(spectra["params"]["charge"][0])).replace("+", "")

    with mgf_output_file as f:

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

def DDA_clean(xl_file: Path,
              mgf_file: Path,
              precursor_delta: float, xlinker_mass: int, ptm_type: str):
    """Clean DDA MGF file.

    Keyword inputs:
    xl_file -- ...
    mgf_file -- ...
    precursor_delta -- ...
    xlinker_type -- ...
    ptm_type -- ...

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.
    """

    # Read the MGF (MS/MS) file.
    mgf_dict_list = list(mgf.read(mgf_file))

    # Importing top_XLs file.
    with open(xl_file, 'r') as f:

        xls = f.read().splitlines()  # Each rows as one element of the list.

    output_file_name = mgf_file.parent / f'{mgf_file.stem}_cleaned.mgf'

    with open(output_file_name, 'w') as f:

        for num_xl, xl in enumerate(xls):

            print(num_xl, xl)

            p1 = ""
            p2 = ""
            fragment_all = []
            mz_light_all = []
            mz_heavy_all = []
            precursor_dict = {}

            (precursor_dict,
             fragment_all,
             mz_light_all,
             mz_heavy_all, p1, p2) = fragment_generator(xl,
                                                        xlinker_mass, ptm_type)

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

    return output_file_name
