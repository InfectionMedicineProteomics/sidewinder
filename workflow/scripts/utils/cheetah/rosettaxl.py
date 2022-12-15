#!/usr/bin/env python

import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import init, Pose
from pyrosetta.rosetta import core, protocols
from pyrosetta.rosetta.protocols.rigid import *

from string import digits
from math import sqrt, pow, exp


## pdb2xl parse the input pdb (according to the partners involved)
## to produce a set of XLs in kojak format according to the distance cut-off.
## inputs: pdb, partners, dist cut-off


def rosetta_eu_dist(pose, a, b):  # CHANGED FROM ALPHA CARBON TO BETA

    a_carbon = 'CA' if pose.residue(a).name1() in 'G' else 'CB'
    b_carbon = 'CA' if pose.residue(b).name1() in 'G' else 'CB'

    dist = sqrt(
        pow((pose.residue(b).atom(b_carbon).xyz().x -
            pose.residue(a).atom(a_carbon).xyz().x), 2)
        + pow((pose.residue(b).atom(b_carbon).xyz().y -
               pose.residue(a).atom(a_carbon).xyz().y), 2)
        + pow((pose.residue(b).atom(b_carbon).xyz().z -
               pose.residue(a).atom(a_carbon).xyz().z), 2)
    )

    return dist


## This function helps to extract each XL peptide from sequence.
## For PEPKPEP--QPEPKQPEP we need to call this function two times to give us
## both peptides. It also helps to define position of K on each peptide which
## is needed in kojak format.
def kojak_generator(sequence, pos):

    xlink = []
    pos_on_pep = 0

    if sequence[pos] == 'K':
        pos_temp1 = 0
        pos_temp2 = 0

        for cnt1 in range(pos-1, -1, -1):
            if (sequence[cnt1] == 'K') or (sequence[cnt1] == 'R'):
                pos_temp1 = cnt1 + 1
                break
            else:
                pos_temp1 = 0

        for cnt2 in range(pos+1, len(sequence), +1):
            if (sequence[cnt2] == 'K') or (sequence[cnt2] == 'R'):
                pos_temp2 = cnt2 + 1
                break
            else:
                pos_temp2 = len(sequence)

        xlink = sequence[pos_temp1:pos_temp2]
        pos_on_pep = pos - pos_temp1 + 1
    else:
        print("ERROR")

    return xlink, pos_on_pep


def rosettaxl(pose, partners, cut_off, out_dir):

    sequence = pose.sequence()
    pdb_info = pose.pdb_info()
    chain_sep = partners.find('_')
    partner1 = partners[0:chain_sep]
    partner2 = partners[chain_sep+1:]

    print("========RosettaXL=========\n")
    print("\nThe protein sequence is: \n", sequence, "\n")

    kojak_list = []
    pep_len = 4

    ALLXL_file = open(out_dir, 'w')

    ## finding and storing all inter_XLs
    for i in range(len(sequence)):

        if sequence[i] == ('K'):
            pep1, K_pos1 = kojak_generator(sequence, i)

            if len(pep1) >= pep_len:

                for j in range(len(sequence)):

                    if sequence[j] == ('K'):
                        pep2, K_pos2 = kojak_generator(sequence, j)

                        if len(pep2) >= pep_len:

                            if pep1 != pep2:

                                if (((pdb_info.chain(i) in partner1) and (pdb_info.chain(j) in partner2)) or ((pdb_info.chain(i) in partner2) and (pdb_info.chain(j) in partner1))):

                                    if (rosetta_eu_dist(pose, i, j) <= cut_off):

                                        XL_kojak_format = "-."+pep1 + \
                                            "("+str(K_pos1)+")--" + \
                                            pep2+"("+str(K_pos2)+").-"
                                        XL_kojak_format_rev = "-."+pep2 + \
                                            "("+str(K_pos2)+")--" + \
                                            pep1+"("+str(K_pos1)+").-"

                                        if (XL_kojak_format not in kojak_list) and (XL_kojak_format_rev not in kojak_list):
                                            kojak_list.append(XL_kojak_format)
                                            print(XL_kojak_format)
                                            ALLXL_file.write(
                                                "%s\n" % XL_kojak_format)
    ALLXL_file.close()
    return
