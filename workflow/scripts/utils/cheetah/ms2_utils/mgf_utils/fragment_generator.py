from pyteomics import mass


def fragments(peptide: str, charge: int):
    """Generate peptide fragments.

    Generate all possible m/z for peptide b and y ions of input charge,
    using the pyteomics mass.fast_mass function to calculate peptide mass.

    Authored by Joel Ströbaek.
    """

    fragment_list = []

    mz_list = []

    ions = ('b', 'y')

    for i in range(1, len(peptide)):

        for ion in ions:

            ion_pep = peptide[:i] if ion in 'b' else peptide[i:]

            pep_mass = mass.fast_mass(ion_pep, ion_type=ion, charge=charge)

            # Adding Carbamidomethylation for Cys residues:
            if "C" in peptide:

                pep_mass += 57.021464 * ion_pep.count("C")

            fragment_list.append(ion_pep)

            mz_list.append(pep_mass)

    return fragment_list, mz_list

def set_ptm(ptm_type, peptide):
    """

    Authored by Joel Ströbaek.
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

def fragment_generator(xl: str, xlinker_mass: int, ptm_type: str):
    """

    Input should be a XL in kojak format.
    example: -.PEPKTIDER(4)--PEPKTIDER(4).-

    Originally authored by Hamed Khakzad, edited by Joel Ströbaek.
    """

    h_mass = 1.008

    fragment_list1 = []

    mz_list1 = []

    fragment_list2 = []

    mz_list2 = []

    #mz_list = []

    fragment_list = []

    comb_frag1 = []

    comb_frag2 = []

    peptide_1, peptide_2 = [peptide.strip('-.(0123456789)')
                            for peptide in xl[0].split('--')]

    p1_k_pos = peptide_1.find('K')

    p2_k_pos = peptide_2.find('K')

    # for filtering spectra based on comparison of mass of XL and the precursor mass
    precursor_dict = {}

    precursor_frag_list = []

    add_to_frag_list_L = []

    add_to_frag_list_H = []

    for precursor_charge in range(3, 9):

        precursor_mz_list_H = []

        precursor_mz_list_L = []

        precursor_charge_letter = ''.join(
            ['+' for s in range(precursor_charge)])

        precursor_frag_list.append(
            peptide_1 + peptide_2 + precursor_charge_letter)

        xl_mass = mass.calculate_mass(sequence=peptide_1) +\
            xlinker_mass +\
            mass.calculate_mass(sequence=peptide_2)

        xl_mass = xl_mass + set_ptm(ptm_type, peptide_1) + \
            set_ptm(ptm_type, peptide_2)

        xl_mass = (xl_mass + (precursor_charge * h_mass))/precursor_charge

        # calculating different isotopic mass/z for the precursor
        precursor_mz_list_L.append(xl_mass)
        precursor_mz_list_L.append(xl_mass + float(1/precursor_charge))
        precursor_mz_list_L.append(xl_mass + float(2/precursor_charge))
        precursor_mz_list_L.append(xl_mass + float(3/precursor_charge))
        precursor_mz_list_L.append(xl_mass + float(4/precursor_charge))

        # calculating M heavy (DSS_D12)
        xl_mass_heavy = xl_mass + ((12*1.006276746)/precursor_charge)
        precursor_mz_list_H.append(xl_mass_heavy)

        precursor_mz_list_H.append(xl_mass_heavy + float(1/precursor_charge))
        precursor_mz_list_H.append(xl_mass_heavy + float(2/precursor_charge))
        precursor_mz_list_H.append(xl_mass_heavy + float(3/precursor_charge))
        precursor_mz_list_H.append(xl_mass_heavy + float(4/precursor_charge))

        # adding the list's to the dictionary
        # 0,1,2,3,4 = Light indices
        # 5,6,7,8,9 = Heavy indices

        precursor_dict[precursor_charge] = (precursor_mz_list_L +
                                           precursor_mz_list_H)

        # here we store light and heavy m/z (each list 6 numbers for different charge values)
        # to add them below to the fragment list
        add_to_frag_list_L.append(xl_mass)

        add_to_frag_list_H.append(xl_mass_heavy)

    # lists to store fragment's and m/z for all charge values
    fragment_list = []

    mz_list_H = []

    mz_list_L = []


    # TODO:
    # - Refactor into function and make function calls instead of below:
    for charge in range(1, 5):

        fragment_list1, mz_list1 = fragments(peptide_1, charge)
        fragment_list2, mz_list2 = fragments(peptide_2, charge)

        charge_letter = ''.join(['+' for s in range(charge)])

        # considering the combination of fragments with full peptide
        ## FRAG1 + XL + PEPTIDE2
        fragment_list_temp1 = list(fragment_list1)

        for num1, frag1 in enumerate(fragment_list_temp1):

            if (peptide_1.find(frag1) <= p1_k_pos) and ((peptide_1.find(frag1) + len(frag1) - 1) >= p1_k_pos):
                fragment_list.append(frag1 + peptide_2 + charge_letter)
                comb_frag1.append(frag1)
                comb_mz = mass.calculate_mass(
                    sequence=frag1) - 18.010565 + xlinker_mass + mass.calculate_mass(sequence=peptide_2)
                #comb_mz = comb_mz + (57.021464*(frag1.count("C")+peptide_2.count("C")))
                comb_mz = comb_mz + \
                    set_ptm(ptm_type, frag1) + set_ptm(ptm_type, peptide_2)

                comb_mz = (comb_mz + (charge * h_mass))/charge

                mz_list_L.append(comb_mz)

                mz_list_H.append(comb_mz + float((12*1.006276746)/charge))

            else:

                fragment_list.append(frag1 + charge_letter)

                mz_list_L.append(mz_list1[num1])

                mz_list_H.append(mz_list1[num1])

        ## FRAG1 + XL + PEPTIDE2
        fragment_list_temp2 = list(fragment_list2)

        for num2, frag2 in enumerate(fragment_list_temp2):

            if (peptide_2.find(frag2) <= p2_k_pos) and (peptide_2.find(frag2) + len(frag2) - 1) >= p2_k_pos:
                fragment_list.append(frag2 + peptide_1 + charge_letter)
                comb_frag2.append(frag2)
                comb_mz = mass.calculate_mass(
                    sequence=frag2) - 18.010565 + xlinker_mass + mass.calculate_mass(sequence=peptide_1)
                #comb_mz = comb_mz + (57.021464*(frag2.count("C")+peptide_1.count("C")))
                comb_mz = comb_mz + \
                    set_ptm(ptm_type, frag2) + set_ptm(ptm_type, peptide_1)

                comb_mz = (comb_mz + (charge * h_mass))/charge

                mz_list_L.append(comb_mz)

                mz_list_H.append(comb_mz + float((12*1.006276746)/charge))

            else:

                fragment_list.append(frag2 + charge_letter)

                mz_list_L.append(mz_list2[num2])

                mz_list_H.append(mz_list2[num2])

        ## FRAG1 + XL + FRAG2
        # here we make all combination of fragments that contain xlinker arm
        # for item1 in comb_frag1:
        #     for item2 in comb_frag2:
        #         if item1+item2 not in fragment_list:
        #             print item1, '+' , item2
        #             fragment_list.append(item1+item2+charge_letter)
        #             comb_mz = mass.calculate_mass(sequence=item1) - 18.010565 + DSS_mass + mass.calculate_mass(sequence=item2) - 18.010565
        #             comb_mz = comb_mz + (57.021464*(item1.count("C")+item2.count("C")))
        #             comb_mz = (comb_mz + (charge * h_mass))/charge
        #             mz_list.append(comb_mz)

    fragment_list.extend(precursor_frag_list)

    mz_list_L.extend(add_to_frag_list_L)

    mz_list_H.extend(add_to_frag_list_H)

    return(precursor_dict,
           fragment_list, mz_list_L, mz_list_H, peptide_1, peptide_2)
