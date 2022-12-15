#!/usr/bin/env python3
"""Simple CLI to generate all possible XLs from an amino acid sequence.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Generate the tryptic peptides separately for the two files and then loop
#   through once and print them (instead of rerunning the kojak generator
#   with every peptide of seq 1); would make the code neater (I think).


from pathlib import Path
from typing import Tuple

import click


def kojak_generator(sequence: str, pos: int) -> Tuple[list, int]:
    """Extract XLs from sequence.

    This function helps to extract each XL peptide from sequence.
    For PEPKPEP--QPEPKQPEP we need to call this function two times to give us
    both peptides. It also helps to define position of K on each peptide which
    is needed in kojak format.

    This function only generates tryptic peptides!

    Originally authored by Hamed Khakzad (Fri Dec  2 14:37:17 2016),
    and edited and commented by Joel Ströbaek (Nov 2022).

    Parameters
    ----------


    Returns
    -------

    """

    xlink = []

    pos_on_pep = 0

    if sequence[pos] in 'K':  # Make sure the correct input position was given.

        pos_temp1 = 0

        pos_temp2 = 0

        # Check upstream for trypsin cleavage site.
        for cnt1 in range(pos-1, -1, -1):

            if (sequence[cnt1] in 'K') or (sequence[cnt1] in 'R'):

                pos_temp1 = cnt1 + 1

                break

            else:

                pos_temp1 = 0

        # Check downstream for trypsin cleavage site.
        for cnt2 in range(pos+1, len(sequence), +1):

            if (sequence[cnt2] in 'K') or (sequence[cnt2] in 'R'):

                pos_temp2 = cnt2 + 1

                break

            else:

                pos_temp2 = len(sequence)

        xlink = sequence[pos_temp1:pos_temp2]

        pos_on_pep = pos - pos_temp1 + 1

    else:

        print("ERROR")

    return xlink, pos_on_pep

def kojak_format(pep_1: str,
                 k_pos_1: int,
                 pep_2: str,
                 k_pos_2: int,
                 swap_order: bool = False) -> str:
    """Generate Kojak formatted string.

    ...

    Parameters
    ----------


    Returns
    -------
    str
    """

    if swap_order:

        return f'-.{pep_2}({str(k_pos_2)})--{pep_1}({str(k_pos_1)}).-'

    else:

        return f'-.{pep_1}({str(k_pos_1)})--{pep_2}({str(k_pos_2)}).-'


@click.command()
@click.version_option(version='1.0')
@click.option(
    '--seq_file1',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help=''
)
@click.option(
    '--seq_file2',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help=''
)
@click.option(
    '--output_dir',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help=''
)
def seq2xl(seq_file1: Path, seq_file2: Path, output_dir: Path) -> Path:
    """Generate all inter-XLs between 2 input aa-seqs.

    Heavily based on code authored by Hamed Khakzad.

    Parameters
    ----------


    Returns
    -------

    """

    pep_len = 4

    kojak_list = []

    seq_1 = ''

    seq_2 = ''

    # Read in sequence text files.
    with open(seq_file1, 'r') as f1, open(seq_file2, 'r') as f2:

        seq_1 = f1.read().strip()

        seq_2 = f2.read().strip()

    xl_file = output_dir / 'all_xls.txt'  # All XLs output file.

    with open(xl_file, 'w') as f:

        # Loop through the seq's and generate all inter-sequence tryptic XLs.
        for i, res_1 in enumerate(seq_1):

            # If the residue is a lysine:
            if res_1 in 'K':

                # Get tryptic peptide + associated kojak format lysine index.
                pep_1, K_pos_1 = kojak_generator(seq_1, i)

                # Only continue if the tryptic peptide is long enough.
                if len(pep_1) >= pep_len:

                    # Check for all potential XLs between the generated
                    # seq 1 tryptic peptide and ones on seq 2.
                    for j, res_2 in enumerate(seq_2):

                        if res_2 in 'K':

                            pep_2, K_pos_2 = kojak_generator(seq_2, j)

                            if len(pep_2) >= pep_len:

                                if pep_1 != pep_2:

                                    # Get kojak formatted string.
                                    kojak_xl = kojak_format(pep_1=pep_1,
                                        k_pos_1=K_pos_1,
                                        pep_2=pep_2,
                                        k_pos_2=K_pos_2)

                                    # Get swapped kojak formatted string.
                                    kojak_xl_swap = kojak_format(pep_1=pep_1,
                                        k_pos_1=K_pos_1,
                                        pep_2=pep_2,
                                        k_pos_2=K_pos_2,
                                        swap_order=True)

                                    # Add unique strings to list and print.
                                    if ((kojak_xl not in kojak_list) and
                                        (kojak_xl_swap not in kojak_list)):

                                        kojak_list.append(kojak_xl)

                                        print(kojak_xl, file=f)

    #return xl_file


if __name__ == '__main__':

    seq2xl()
