from string import digits
from math import pow, exp

from Bio.PDB import CaPPBuilder, Selection

## xlvalidation parse the input pdb (according to the partners involved)
## to obtain number of XLs it fulfills.
## inputs: pdb, list of XLs, partners, dist cut-off

def eu_dist(structure, a, b):

    chain1 = Selection.unfold_entities(structure, 'R')[a].get_parent().id
    chain2 = Selection.unfold_entities(structure, 'R')[b].get_parent().id

    atom1 = structure[0][chain1][a]['CA']
    atom2 = structure[0][chain2][b]['CA']

    return atom2 - atom1


def xlvalidation(structure, top_XL_file, cut_off):

    ## 1. Reading and storing the pdb in a pose
    #parser = PDBParser()
    #structure = parser.get_structure('INPS', './decoy1.pdb')
    #model = structure[0]

    sequence_list_A = []
    sequence_list_B = []
    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure[0]['A']):
        sequence_list_A.append(str(pp.get_sequence()))
    sequence_A = ''.join(sequence_list_A)
    print ("\nPartner 1 sequence is: \n", sequence_A)

    for pp in ppb.build_peptides(structure[0]['B']):
        sequence_list_B.append(str(pp.get_sequence()))

    sequence_B = ''.join(sequence_list_B)
    print ("\nPartner 2 sequence is: \n", sequence_B)

    sequence_combined = sequence_A + sequence_B
    print ("\nCombined sequence is: \n", sequence_combined)


    ## 2. Reading the xl_file
    with open(top_XL_file,'r') as xlfile:
        top_XL = xlfile.read().splitlines() # each rows as one element of the list

    output_xl_number = 0
    normal_dist_score = 0.0
    good_XL_list = []
    for num_xl,xl in enumerate(top_XL):
        print (num_xl+1, xl)

        K_pos_P1 = 0
        K_pos_P2 = 0
        eulidean_dist = 10000.0

        xl_trans = str.maketrans('', '', digits)
        xl_without_digit = xl.translate(xl_trans)

        peptide1 = xl_without_digit.split('--')[0].replace("-","").replace(".","")[:-2]
        peptide2 = xl_without_digit.split('--')[1].replace("-","")[:-3]
        K_pos_P1 = peptide1.find('K') + 1
        K_pos_P2 = peptide2.find('K') + 1

        multiple_occ_list1 = [x for x in range(len(sequence_combined)) if sequence_combined.find(peptide1, x) == x]
        multiple_occ_list2 = [y for y in range(len(sequence_combined)) if sequence_combined.find(peptide2, y) == y]
        ## finding minimum distance if multiple occurance happened
        seq_pos_p1_k = 0
        seq_pos_p2_k = 0
        tmp_dist = 10000.0
        try:
            for pos1 in multiple_occ_list1:
                for pos2 in multiple_occ_list2:
                    if tmp_dist > eu_dist(structure, pos1+K_pos_P1+1, pos2+K_pos_P2+1):
                        tmp_dist = eu_dist(structure, pos1+K_pos_P1+1, pos2+K_pos_P2+1)
                        seq_pos_p1_k = pos1+K_pos_P1
                        seq_pos_p2_k = pos2+K_pos_P2

            print ("aa positions on the sequences: ", seq_pos_p1_k, seq_pos_p2_k)

            if ((peptide1 in sequence_combined) and (peptide2 in sequence_combined)):
                eulidean_dist = eu_dist(structure, seq_pos_p1_k+1, seq_pos_p2_k+1)

                print ("Euclidean distance is:  ", eulidean_dist, "\n")

            else:
                print ("The XL is not found on the protein sequence. Check each peptide to be valid!\n")
                # out_sc_file.write("XL is not found\n")


            if eulidean_dist <= cut_off:
                good_XL_list.append(xl)
                output_xl_number += 1
                normal_dist_score += (10 * (1/(exp ( pow( (eulidean_dist-17.5),2 )/30 )))) #Normal Distribution Formula Ae^-(x-M)/2q^2
                # out_sc_file.write(xl+", ")
                print ("found a valid XL!", eulidean_dist, "\n")
        except:
            pass

    # out_sc_file.write(str(output_xl_number)+", "+str(normal_dist_score))
    print ("Score:", output_xl_number, normal_dist_score)
    return output_xl_number, normal_dist_score, good_XL_list
