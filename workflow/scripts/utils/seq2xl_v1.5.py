#!/usr/bin/env python3
"""Basic CLI to generate all possible XLs from an amino acid sequence.
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
from Bio.PDB import PDBConstructWarning


warnings.simplefilter('ignore', BiopythonParserWarning)

warnings.simplefilter('ignore', PDBConstructWarning)


def xl_generator(sequence: str,
                 pos: int,
                 xl_residue: str = 'K',
                 len_cutoff: int = 4) -> str:
    """Extracts a potential XL peptide from a given sequence position.

    This function generates a candidate XL peptide from a protein sequence,
    considering a specific position and the provided XL-linkable residue type.
    It assumes tryptic digestion and only extracts peptides between upstream
    and downstream Lysine (K) or Arginine (R) residues.
    For PEPKPEP--QPEPKQPEP we need to call this function two times to give us
    both peptides.

    Originally authored by Hamed Khakzad (Fri Dec  2 14:37:17 2016),
    and edited and commented by Joel Ströbaek (Nov 2022).

    Args:
        sequence (str): The protein sequence.
        pos (int): The zero-based index of the potential XL-linkable residue.
        xl_residue (str, optional): The amino acid type considered as the
            XL-linkable residue. Defaults to 'K' (Lysine).
        len_cutoff (int, optional): The minimum allowed length for a valid
            XL peptide. Defaults to 4.

    Returns:
        str: The extracted XL peptide sequence if it meets the length cutoff,
            otherwise None.

    Raises:
        ValueError: If the provided position is outside the sequence length
            or the specified residue at that position is not the XL-linkable type.
    """
    #TODO
    #   - Make (or implement) a proper digestor

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
    """Finds the zero-based index of the first XL-linkable residue in a sequence.

    This function searches for the first occurrence of the specified XL-linkable
    residue within a protein sequence and returns its index. The indexing can be
    adjusted based on the `kojak_indexing` flag.

    Args:
        seq (str): The protein sequence.
        xl_residue (str, optional): The amino acid type considered as the
            XL-linkable residue. Defaults to 'K' (Lysine).
        kojak_indexing (bool, optional): Flag to adjust the indexing by 1
            (True) for compatibility with Kojak output, or use zero-based
            indexing (False). Defaults to True.

    Returns:
        int: The zero-based (or adjusted) index of the first XL-linkable residue
            in the sequence, or -1 if not found.
    """
    kojak_index = 1 if kojak_indexing else 0

    return seq.find(xl_residue) + kojak_index

def kojak_format(pep_1: str,
                 k_pos_1: int,
                 pep_2: str,
                 k_pos_2: int,
                 swap_order: bool = False) -> str:
    """Formats an XL interaction string in Kojak format.

    This function generates a string representing an XL crosslink between two
    peptides in the format used by the Kojak software. It allows specifying
    the peptide sequences, their respective XL-linkable residue positions, and
    an optional order swap for the peptides.

    Args:
        pep_1 (str): The sequence of the first peptide.
        k_pos_1 (int): The zero-based index of the XL-linkable residue in pep_1.
        pep_2 (str): The sequence of the second peptide.
        k_pos_2 (int): The zero-based index of the XL-linkable residue in pep_2.
        swap_order (bool, optional): Flag to swap the order of peptides in the
            output string (True), or keep the original order (False). Defaults
            to False.

    Returns:
        str: The formatted Kojak string representing the XL interaction.
    """
    if swap_order:

        return f'-.{pep_2}({str(k_pos_2)})--{pep_1}({str(k_pos_1)}).-'

    else:

        return f'-.{pep_1}({str(k_pos_1)})--{pep_2}({str(k_pos_2)}).-'

def protein_to_rec_list(pdb_path: Path,
                        file_type: str = 'pdb-atom') -> List[SeqIO.SeqRecord]:
    """Parses a PDB file and returns a list of SeqRecord objects.

    This function reads a Protein Data Bank (PDB) file and utilizes Biopython's
    SeqIO module to parse the protein sequences into SeqRecord objects.

    Args:
        pdb_path (Path): Path to the PDB file.
        file_type (str, optional): Format specification for SeqIO parsing.
            Defaults to 'pdb-atom'.

    Returns:
        List[SeqIO.SeqRecord]: A list containing SeqRecord objects representing
            the protein sequences from the PDB file.

    Raises:
        FileNotFoundError: If the specified PDB file cannot be found.
    """
    chain_recs = []

    with open(pdb_path, 'r+') as f:

        for rec in SeqIO.parse(f, file_type):

            chain_recs.append(rec)

    return chain_recs

def rec_to_xls(rec_list: List[SeqIO.SeqRecord],
               linker_res: str = 'K',
               max_peptide_length: int = 4) -> List[tuple]:
    """Extracts potential XL-linkable peptides from a list of SeqRecord objects.

    This function iterates through a list of SeqRecord objects (representing
    protein sequences) and identifies potential XL peptides based on the
    provided linker residue type. It considers tryptic digestion and extracts
    peptides between upstream and downstream cleavage sites (Lysine or Arginine).
    Only peptides exceeding a minimum length cutoff are retained as candidate
    XL peptides.

    Args:
        rec_list (List[SeqIO.SeqRecord]): A list containing SeqRecord objects.
        linker_res (str, optional): The amino acid type considered as the
            XL-linkable residue. Defaults to 'K' (Lysine).
        max_peptide_length (int, optional): The minimum allowed length for a
            valid XL peptide. Defaults to 4.

    Returns:
        List[tuple]: A list of tuples where each tuple represents a unique
            XL peptide candidate, containing the peptide sequence and the
            zero-based index of the XL-linkable residue within that sequence.
    """
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
    """Generates possible inter-XL interactions between two PDB files.

   This function takes two PDB files as input and generates potential
   inter-XL interactions between their protein sequences. It considers
   tryptic digestion and filters peptides based on a minimum length cutoff.
   Optional validation of XLs from the second PDB file can be performed
   using a FASTA file for a superset structure. Results are written to
   an output file in Kojak format.

   Heavily based on code authored by Hamed Khakzad.

   Args:
       pdb_file1 (Path): Path to the first PDB file.
       pdb_file2 (Path): Path to the second PDB file.
       output_dir (Path): Path to the output directory.
       max_peptide_length (int): Minimum allowed length for XL peptides.
       validation_fasta (Path): Path to a FASTA file for validation (optional).

   Raises:
       ValueError: If validation_fasta is provided but not a valid FASTA file.
           FileNotFoundError: If input PDB files cannot be found.
           IOError: If output file cannot be written.
   """
    pep_len = max_peptide_length

    kojak_xls = []

    # Read in PDB to SeqRecord:
    rec_1 = protein_to_rec_list(pdb_file1)

    rec_2 = protein_to_rec_list(pdb_file2)

    # Extract XLs from SeqRecords:
    xls_1 = rec_to_xls(rec_1, max_peptide_length=pep_len)

    xls_2 = rec_to_xls(rec_2, max_peptide_length=pep_len)

    # All XLs output file.
    xl_file = output_dir / f'{pdb_file2.stem}_+_{pdb_file1.stem}.xls'

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
