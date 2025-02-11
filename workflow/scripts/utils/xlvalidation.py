#!/usr/bin/env python3

from string import digits
from math import pow, exp
from typing import List, Tuple

from Bio.PDB import CaPPBuilder, Selection, Structure


def eu_dist(structure: Structure, i: int, j: int) -> float:
    """Calculates the Euclidean distance between two alpha-carbon (CA) atoms within a protein structure.

    This function assumes a PDB structure containing protein chains with a 'CA' atom for each residue.
    It retrieves the CA atoms for the residues specified by indices `i` and `j`, regardless of their chains,
    and calculates the Euclidean distance between them.

    Args:
        structure (Structure): A Bio.PDB Structure object containing the protein structure.
        i (int): Index of the first residue for distance calculation.
        j (int): Index of the second residue for distance calculation.

    Returns:
        float: The Euclidean distance between the CA atoms of residues i and j

    Raises:
        ValueError: If the structure or residues do not meet expected format or content.
    """
    chain1 = Selection.unfold_entities(structure, 'R')[i].get_parent().id

    chain2 = Selection.unfold_entities(structure, 'R')[j].get_parent().id

    atom1 = structure[0][chain1][i]['CA']

    atom2 = structure[0][chain2][j]['CA']

    return atom2 - atom1

def xlvalidation(structure: Structure,
                 top_xls: List[str],
                 cut_off: float) -> Tuple[int, float, List[str]]:
    """
    Validates cross-links (XLs) within a protein structure and calculates a score based on Euclidean distances.

    This function analyzes a protein structure (PDB format) and a list of top cross-links (XLs)
    to determine how many XLs are satisfied based on a distance cutoff. It also calculates a score
    for each valid XL using a normal distribution formula and the Euclidean distance between the
    cross-linked residues' Cα atoms.

    Args:
        structure (Structure): A Bio.PDB Structure object representing the protein structure.
        top_xls (List[str]): List of strings containing the top cross-link definitions.
            Each string should follow a specific format (refer to documentation).
        cut_off (float): Distance cutoff for considering a cross-link as satisfied.

    Returns:
        Tuple[int, float, List[str]]:
            - int: Number of satisfied cross-links.
            - float: Total score based on Euclidean distances for valid cross-links.
            - List[str]: List of valid cross-link definitions from the input.

    Raises:
        ValueError: If the structure, XL definitions, or distance cutoff are invalid.
    """
    sequence_list_A = []

    sequence_list_B = []

    ppb = CaPPBuilder()

    for pp in ppb.build_peptides(structure[0]['A']):

        sequence_list_A.append(str(pp.get_sequence()))

    sequence_A = ''.join(sequence_list_A)

    for pp in ppb.build_peptides(structure[0]['B']):

        sequence_list_B.append(str(pp.get_sequence()))

    sequence_B = ''.join(sequence_list_B)

    sequence_combined = sequence_A + sequence_B

    output_xl_number = 0

    normal_dist_score = 0.0

    good_XL_list = []

    for num_xl, xl in enumerate(top_xls):

        K_pos_P1 = 0

        K_pos_P2 = 0

        euclidean_dist = 10000.0

        xl_trans = str.maketrans('', '', digits)

        xl_without_digit = xl.translate(xl_trans)

        peptide1 = xl_without_digit.split('--')[0].replace("-","").replace(".","")[:-2]

        peptide2 = xl_without_digit.split('--')[1].replace("-","")[:-3]

        K_pos_P1 = peptide1.find('K') + 1

        K_pos_P2 = peptide2.find('K') + 1

        multiple_occ_list1 = [x for x in range(len(sequence_combined))
                              if sequence_combined.find(peptide1, x) == x]

        multiple_occ_list2 = [y for y in range(len(sequence_combined))
                              if sequence_combined.find(peptide2, y) == y]

        ## Finding minimum distance if multiple occurance happened.
        seq_pos_p1_k = 0

        seq_pos_p2_k = 0

        tmp_dist = 10000.0

        try:

            for pos1 in multiple_occ_list1:

                for pos2 in multiple_occ_list2:

                    if tmp_dist > eu_dist(structure, pos1+K_pos_P1+1,
                                          pos2+K_pos_P2+1):

                        tmp_dist = eu_dist(structure, pos1+K_pos_P1+1,
                                           pos2+K_pos_P2+1)

                        seq_pos_p1_k = pos1+K_pos_P1

                        seq_pos_p2_k = pos2+K_pos_P2

            if ((peptide1 in sequence_combined) and
                (peptide2 in sequence_combined)):

                euclidean_dist = eu_dist(structure, seq_pos_p1_k+1,
                                        seq_pos_p2_k+1)

            if euclidean_dist <= cut_off:

                good_XL_list.append(xl)

                output_xl_number += 1

                d_opti = cut_off*.55  # Hypothetical distance optimum.

                # Normal distribution formula: Ae^-(x-M)/2q^2
                normal_dist_score += (10*(1/(exp(pow((euclidean_dist-d_opti),
                                                       2)/cut_off))))

        except:

            pass

    return (output_xl_number, normal_dist_score, good_XL_list)
