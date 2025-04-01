#!/usr/bin/env python3
from typing import Dict

from pandas import DataFrame, Series

from . import crosslinking


def add_xl_info(x: Series,
                ab_chain_dict: Dict,
                ag_section_dict: Dict,
                is_dimer: bool = False) -> Series:

    hits = []

    xl = x['XL']

    ab_pep = xl.strip('.-').split('--')[0].split('(')[0]

    hc = x['heavy_chain'][:-2]

    hc_len = len(ab_chain_dict[hc]['seq_record'])

    lc = x['light_chain'][:-2]

    lc_len = len(ab_chain_dict[lc]['seq_record'])

    ag_start = hc_len + lc_len

    # source = hc if ab_chain_dict[hc]['seq_record'].find(ab_pep) >= 0 else lc

    source_list = [chain for chain in [hc, lc]
                   if ab_pep in str(ab_chain_dict[chain]['seq_record'])]

    x['XL_source'] = ','.join(source_list)

    for source in source_list:

        ab_seq = str(ab_chain_dict[source]['seq_record'])

        ag_seq = ag_section_dict[x['antigen']]

        hits.extend(crosslinking.xl_finder(xl, p1_seq=ab_seq, p2_seq=ag_seq))

        ab_start = 0 if source == lc else lc_len

    hit = [(y[0] + ab_start, y[1] + (ag_start), y[1] + (ag_start) + len(ag_seq))
           if is_dimer
           else (y[0] + ab_start, y[1] + (ag_start))
           for y in hits]

    x['pyMOL_XL_pos'] = hit

    return x

def peptide_hits(x: Series) -> Series:

    xl = x['XL']

    p1, p2 = [p.split('(')[0]
              for p in xl.strip('-.').split('--')]

    ions = [ion.rstrip('+') for ion in x['covered_Frags'].split(',')]

    ions = list(set(ions))

    hits = [0, 0]

    for ion in ions:

        if (p1.startswith(ion)) or (p1.endswith(ion)):

            hits[0] += 1

        elif (p2.startswith(ion)) or (p2.endswith(ion)):

            hits[1] += 1

        else:

            continue

    good_assignment = False

    if (hits[0] >= 3) & (hits[1] >= 3):

        good_assignment = True

    x['hits'] = hits

    x['good_assignment'] = good_assignment

    return x


def score(x: Series) -> Series:

    ions = [ion.rstrip('+') for ion in x['covered_Frags'].split(',')]

    xl = x['XL']

    p1, p2 = [p.split('(')[0]
              for p in xl.strip('-.').split('--')]

    p1_score = 0

    p2_score = 0

    xl_score = 0

    components = 0

    for ion in ions:

        a = p1.find(ion)

        b = p2.find(ion)

        if (a >= 0) or (b >= 0):

            if (p1[a:] == ion) | (a == 0):

                p1_score += len(ion) / len(p1)

                components += 1

            if (p2[b:] == ion) | (b == 0):

                p2_score += len(ion) / len(p2)

                components += 1

        else:

            xl_score += len(ion) / (len(p1) + len(p2))

            components += 1

    score = (components / 3) * (p1_score + p2_score + xl_score)

    x['spectra_score'] = score

    return x

def min_max_norm_col(dataframe: DataFrame, col_selection: str) -> DataFrame:

    df = dataframe.copy()

    df[col_selection] = df[col_selection] / df[col_selection].max()

    return df
