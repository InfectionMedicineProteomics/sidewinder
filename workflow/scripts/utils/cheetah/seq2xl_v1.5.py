#!/usr/bin/env python3
"""Simple CLI to generate all possible XLs from an amino acid sequence.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


#TODO
# - Adjust for non-tryptic digest
# - Make validation available for both input proteins.


import warnings
from pathlib import Path
from typing import List

import click
from Bio import BiopythonParserWarning, SeqIO

warnings.simplefilter('ignore', BiopythonParserWarning)


def xl_generator(sequence: str,
                 pos: int,
                 xl_residue: str = 'K',
                 len_cutoff: int = 4) -> str:
    """Extract XLs from sequence.

    This function helps to extract each XL peptide from sequence.
    For PEPKPEP--QPEPKQPEP we need to call this function two times to give us
    both peptides.

    This function only generates tryptic peptides!

    Originally authored by Hamed Khakzad (Fri Dec  2 14:37:17 2016),
    and edited and commented by Joel Ströbaek (Nov 2022).

    Parameters
    ----------


    Returns
    -------

    """

    #TODO
    #   - Make a proper digestor

    xlink = []

    # Make sure the correct input position was given.
    if sequence[pos] in xl_residue:

        pos_temp1 = 0

        pos_temp2 = 0

        # Check upstream for trypsin cleavage site.
        for cnt1 in range(pos - 1, -1, -1):

            if (sequence[cnt1] in 'K') or (sequence[cnt1] in 'R'):

                pos_temp1 = cnt1 + 1

                break

            else:

                pos_temp1 = 0

        # Check downstream for trypsin cleavage site.
        for cnt2 in range(pos + 1, len(sequence), +1):

            if (sequence[cnt2] in 'K') or (sequence[cnt2] in 'R'):

                pos_temp2 = cnt2 + 1

                break

            else:

                pos_temp2 = len(sequence)

        xlink = sequence[pos_temp1:pos_temp2]

    else:

        print("ERROR")

    if len(xlink) >= len_cutoff:

        return xlink

    else:

        return None

def xl_residue_index(seq: str,
                     xl_residue: str = 'K', kojak_indexing: bool = True):

    kojak_index = 1 if kojak_indexing else 0

    return seq.find(xl_residue) + kojak_index

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

def protein_to_rec_list(pdb_path: Path, file_type: str = 'pdb-atom'):

    chain_recs = []

    with open(pdb_path, 'r+') as f:

        for rec in SeqIO.parse(f, file_type):

            chain_recs.append(rec)

    return chain_recs

def rec_to_xls(rec_list,
               linker_res: str = 'K',
               max_peptide_length: int = 4) -> List[tuple]:

    prev_chains = []

    xls = []

    for rec in rec_list:

        seq = str(rec.seq)

        if seq in prev_chains:

            continue

        else:

            prev_chains.append(seq)

            for i, r in enumerate(seq):

                if r in linker_res:

                    xl = xl_generator(seq, i, linker_res, max_peptide_length)

                    if xl:

                        xl_tuple = (xl, xl_residue_index(xl))

                        if xl_tuple not in xls:

                            xls.append(xl_tuple)

    return xls


@click.command()
@click.version_option(version='1.5')
@click.option('--pdb_file1',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--pdb_file2',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--output_dir',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--validation_fasta',
              required=False,
              type=click.Path(exists=True, path_type=Path),
              default=None,
              help=str(f'Used to validate the XLs assigned from pdb_file2, '
                       f'if the associated PDB contains a subset structure '
                       f'from the validation_pdb.'))
@click.option('--max_peptide_length',
              required=False,
              type=int,
              default=4, help='Ignore [tryptic] peptides shorter than this.')
def seq2xl(pdb_file1: Path,
           pdb_file2: Path,
           output_dir: Path,
           max_peptide_length: int, validation_fasta: Path):
    """Generate all inter-XLs between 2 input aa-seqs.

    Heavily based on code authored by Hamed Khakzad. Only considers
    tryptic digests in current iteration.

    Parameters
    ----------


    Returns
    -------

    """

    pep_len = max_peptide_length

    kojak_xls = []

    # Read in PDB to SeqRecord:
    rec_1 = protein_to_rec_list(pdb_file1)

    rec_2 = protein_to_rec_list(pdb_file2)

    # Extract XLs from SeqRecords:
    xls_1 = rec_to_xls(rec_1, max_peptide_length=pep_len)

    xls_2 = rec_to_xls(rec_2, max_peptide_length=pep_len)

    xl_file = output_dir / 'all_xls.txt'  # All XLs output file.

    # Validate the XLs assigned to xls_2 based on XLs generated for a
    # superset structure of the protein. Useful if not the entire protein
    # is being targeted (to make sure there are no in silico artifacts).
    if validation_fasta:

        val_rec = protein_to_rec_list(validation_fasta, 'fasta')

        val_xls = rec_to_xls(val_rec, max_peptide_length=pep_len)

        val_xl_ref = [x[0] for x in val_xls]

        for xl, idx in xls_2:

            if xl not in val_xl_ref:

                xls_2.remove((xl, idx))

    for pep_1, idx_1 in xls_1:

        for pep_2, idx_2 in xls_2:

            kojak_xl = kojak_format(pep_1, idx_1, pep_2, idx_2)

            kojak_xl_swap = kojak_format(pep_1,
                                         idx_1, pep_2, idx_2, swap_order=True)

            if (kojak_xl not in kojak_xls) and (kojak_xl_swap not in kojak_xls):

                kojak_xls.append(kojak_xl)

    with open(xl_file, '+w') as out_file:

        for xl in kojak_xls:

            out_file.write(f'{xl}\n')


if __name__ == '__main__':

    seq2xl()
