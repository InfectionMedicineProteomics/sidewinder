#!/usr/bin/env python3

# TODO:
# - Look into comment:
#       "# Calculating M heavy (DSS_D12)", is this hardcoding DSS despite prev.
#       script allowing for setting of the cross-linker?


from typing import Tuple, List
from pyteomics import mass


def fragments(peptide: str, charge: int) -> Tuple[List[str], List[float]]:
    """Generate peptide fragments.

    This function generates all possible b and y ion fragments for a given
    peptide sequence and charge state. It utilizes the
    `pyteomics.mass.fast_mass` function to calculate the mass of each fragment,
    considering the mass of cysteine modifications (57.021464 Da).

    Args:
        peptide (str): The amino acid sequence of the peptide.
        charge (int): The charge state of the peptide.

    Returns:
        Tuple[List[str], List[float]]: A tuple containing two lists:
            - A list of fragment ion names (in the order of b1, y1, b2, ...,
            yN).
            - A list of corresponding m/z values for the fragments.
    """
    b_frag = [peptide[:i] for i in range(1, len(peptide))]

    b_ions = [mass.fast_mass(frag,
              ion_type='b', charge=charge) + 57.021464 * frag.count("C")
              for frag in b_frag]

    y_frag = [peptide[i:] for i in range(1, len(peptide))]

    y_ions = [mass.fast_mass(frag,
              ion_type='y', charge=charge) + 57.021464 * frag.count("C")
              for frag in y_frag]

    return ([x for tuples in zip(b_frag, y_frag) for x in tuples],
            [x for tuples in zip(b_ions, y_ions) for x in tuples])

def calc_ptm_mass(ptm_type: str, peptide: str) -> float:
    """Calculates the mass shift due to a PTM applied to a peptide.

    This function takes a PTM type (as a string) and a peptide sequence and
    returns the total mass difference caused by the PTMs on that peptide. It
    uses a predefined dictionary to map PTM types to their mass shifts and
    amino acid targets.

    Currently, this function only supports unmodified peptides (ptm_type="1").

    Args:
        ptm_type (str): The type of PTM applied (e.g., "1" for unmodified).
        peptide (str): The amino acid sequence of the peptide.

    Returns:
        float: The total mass shift caused by PTMs on the peptide (0 for
            unmodified).
    """
    ptm_dict = {"1": (57.021464, "C"),
                "2": (-79.966, "Y"),
                "3": (-79.966, "T"),
                "4": (-79.966, "S"),
                "5": (42.011, "K"),
                "6": (14.016, "K"),
                "7": (28.031, "K"),
                "8": (42.047, "K"),
                "9": (57.021464, "R"),
                "10": (143.155, "R"),
                "11": (15.995, "M"),
                "12": (3.995, "W")}

    ptm_mass = ptm_dict[ptm_type][0]

    ptm_count = peptide.count(ptm_dict[ptm_type][1])

    return 0 + (ptm_mass * ptm_count)

def fragment_generator(xl: str, xlinker_mass: int, xlinker: int, ptm_type: str) -> Tuple:
    """Generates theoretical fragment ions for a cross-linked peptide.

    This function takes a cross-linked peptide in Kojak format (e.g.,
    ".--PEPTIDE(1)--PEPTIDE(2)."), the cross-linker mass, and the PTM type
    (currently supports unmodified peptides only), and returns a tuple
    containing various information for analysis:

    - A dictionary mapping precursor charge states to their theoretical m/z values (light and heavy isotopes).
    - A list of fragment ion sequences (considering all combinations of peptide fragments and the cross-linker).
    - Two lists containing m/z values for the light and heavy isotope forms of each fragment ion.
    - The amino acid sequences of the two peptides linked by the cross-linker.

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.

    Args:
        xl (str): The cross-linked peptide sequence in Kojak format.
        xlinker_mass (int): The mass of the cross-linker used.
        ptm_type (str): The type of PTM applied (currently supports unmodified,
            "1").

    Returns:
        Tuple: A tuple containing the following elements:
            - precursor_dict (dict): Maps precursor charge states to their m/z values (light and heavy).
            - fragment_list (List[str]): List of fragment ion sequences.
            - mz_list_L (List[float]): m/z values for light isotope fragment ions.
            - mz_list_H (List[float]): m/z values for heavy isotope fragment ions.
            - peptide_1 (str): Amino acid sequence of the first peptide.
            - peptide_2 (str): Amino acid sequence of the second peptide.
    """
    h_mass = 1.008

    xlinker_deut = {1: 12*1.006276746,  # DSS
                    2: 6*1.006276746,  # DSG
                    3: 12*1.006276746}  # EGS

    mz_list_1 = []

    mz_list_2 = []

    fragment_list = []

    peptide_1, peptide_2 = [peptide.strip('-.(0123456789)')
                            for peptide in xl.split('--')]

    p1_k_pos = peptide_1.find('K')

    p2_k_pos = peptide_2.find('K')

    p1_mass = mass.calculate_mass(sequence=peptide_1)

    p2_mass = mass.calculate_mass(sequence=peptide_2)

    p1_ptms = calc_ptm_mass(ptm_type, peptide_1)

    p2_ptms = calc_ptm_mass(ptm_type, peptide_2)

    # Filtering spectra based on mass comparison for XL and the precursor mass
    precursor_dict = {}

    precursor_frag_list = []

    add_to_frag_list_L = []

    add_to_frag_list_H = []

    for precursor_charge in range(3, 9):

        precursor_charge_letter = '+' * precursor_charge

        precursor_frag_list.append(peptide_1
                                   + peptide_2 + precursor_charge_letter)

        tmp_mass = p1_mass + xlinker_mass + p2_mass + p1_ptms + p2_ptms

        xl_mass = (tmp_mass + (precursor_charge * h_mass))/precursor_charge

        # Calculating different isotopic m/z for the precursor:
        precursor_mz_list_L = [xl_mass + float(i / precursor_charge)
                               for i in range(5)]

        # Calculating M heavy
        xl_mass_heavy = xl_mass + ((xlinker_deut[xlinker])/precursor_charge)

        precursor_mz_list_H = [xl_mass_heavy + float(i / precursor_charge)
                               for i in range(5)]

        # Adding the list's to the dictionary:
        # 0,1,2,3,4 = Light indices
        # 5,6,7,8,9 = Heavy indices
        precursor_dict[precursor_charge] = (precursor_mz_list_L
                                            + precursor_mz_list_H)

        # Store light and heavy m/z (each list 6 numbers for different
        # charge values) to add them below to the fragment list.
        add_to_frag_list_L.append(xl_mass)

        add_to_frag_list_H.append(xl_mass_heavy)

    # Lists to store fragment's and m/z for all charge values.
    fragment_list = []

    mz_list_H = []

    mz_list_L = []


    # TODO:
    # - Refactor into function and make function calls instead of below:
    for charge in range(1, 5):

        fragment_1_list, mz_list_1 = fragments(peptide_1, charge)
        fragment_2_list, mz_list_2 = fragments(peptide_2, charge)

        charge_letter = '+' * charge

        # Considering the combination of fragments with full peptide.
        # Peptide 1 fragment + XL + peptide 2:
        for num1, frag1 in enumerate(fragment_1_list):

            if ((peptide_1.find(frag1) <= p1_k_pos)
                and ((peptide_1.find(frag1) + len(frag1) - 1) >= p1_k_pos)):

                f1_mass = mass.calculate_mass(sequence=frag1)

                f1_ptms = calc_ptm_mass(ptm_type, frag1)

                fragment_list.append(frag1 + peptide_2 + charge_letter)

                #comb_frag1.append(frag1)

                tmp_comb_mz = (f1_mass
                               - 18.010565
                               + xlinker_mass + p2_mass + f1_ptms + p2_ptms)
                #comb_mz = comb_mz + (57.021464
                #                     * (frag1.count("C")
                #                        + peptide_2.count("C")))

                comb_mz = (tmp_comb_mz + (charge * h_mass))/charge

                mz_list_L.append(comb_mz)

                mz_list_H.append(comb_mz + float((12 * 1.006276746)/charge))

            else:

                fragment_list.append(frag1 + charge_letter)

                mz_list_L.append(mz_list_1[num1])

                mz_list_H.append(mz_list_1[num1])

        # Peptide 2 fragment + XL + peptide 1:
        for num2, frag2 in enumerate(fragment_2_list):

            if ((peptide_2.find(frag2) <= p2_k_pos)
                and (peptide_2.find(frag2) + len(frag2) - 1) >= p2_k_pos):

                f2_mass = mass.calculate_mass(sequence=frag2)

                f2_ptms = calc_ptm_mass(ptm_type, frag2)

                fragment_list.append(frag2 + peptide_1 + charge_letter)

                #comb_frag2.append(frag2)

                tmp_comb_mz = (f2_mass
                               - 18.010565
                               + xlinker_mass + p1_mass + f2_ptms + p1_ptms)
                #comb_mz = comb_mz + (57.021464
                #                     * (frag2.count("C")
                #                        + peptide_1.count("C")))

                comb_mz = (tmp_comb_mz + (charge * h_mass))/charge

                mz_list_L.append(comb_mz)

                mz_list_H.append(comb_mz + float((12 * 1.006276746)/charge))

            else:

                fragment_list.append(frag2 + charge_letter)

                mz_list_L.append(mz_list_2[num2])

                mz_list_H.append(mz_list_2[num2])

        # Peptide 1 fragment + XL + peptide 2 fragment:
        # Here we make all combination of fragments that contain the
        # cross-linker arm.
        # for item1 in comb_frag1:
        #     for item2 in comb_frag2:
        #         if item1 + item2 not in fragment_list:
        #             i1_mass = mass.calculate_mass(sequence=item1)
        #             i2_mass = mass.calculate_mass(sequence=item2)
        #             i1_c = item1.count("C")
        #             i2_c = item2.count("C")
        #             fragment_list.append(item1 + item2 + charge_letter)
        #             tmp_comb_mz = (i1_mass
        #                            + xlinker_mass
        #                            - (18.010565 * 2)
        #                            + i2_mass) + (57.021464 * (i1_c + i2_c))
        #             comb_mz = (tmp_comb_mz + (charge * h_mass))/charge
        #             mz_list.append(comb_mz)

    fragment_list.extend(precursor_frag_list)

    mz_list_L.extend(add_to_frag_list_L)

    mz_list_H.extend(add_to_frag_list_H)

    return(precursor_dict,
           fragment_list, mz_list_L, mz_list_H, peptide_1, peptide_2)
