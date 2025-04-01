#!/usr/bin/env python3
from typing import List, Tuple

from Bio.PDB import Selection, Structure


def eu_dist(structure: Structure, a: int, b: int) -> float:

    chain_1 = Selection.unfold_entities(structure, 'R')[a].get_parent().id

    chain_2 = Selection.unfold_entities(structure, 'R')[b].get_parent().id

    atom_1 = structure[0][chain_1][a]['CA']

    atom_2 = structure[0][chain_2][b]['CA']

    return atom_2 - atom_1

def xl_finder(xl, p1_seq, p2_seq) -> List[Tuple]:

    p1, p2 = xl.strip('-.').split('--')

    p1, k1 = p1.strip(')').split('(')

    p2, k2 = p2.strip(')').split('(')

    p1_pos_list = [x for x in range(len(p1_seq))
                   if p1_seq.find(p1, x) == x]

    p2_pos_list = [x for x in range(len(p2_seq))
                   if p2_seq.find(p2, x) == x]

    pos_list = [(x + int(k1), y + int(k2))
                for x in p1_pos_list
                for y in p2_pos_list]

    return pos_list
