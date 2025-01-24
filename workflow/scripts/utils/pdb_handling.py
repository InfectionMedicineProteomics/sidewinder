#!/usr/bin/env python3
"""Extract AA sequence, and clean/format input PDB files.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


# TODO:
# - Potentially change out CLI pdb_tools for pdbfixer python package.
# - Make sure invalid amino acids are caught!


from pathlib import Path
from subprocess import Popen, PIPE
from typing import List

from Bio.PDB import PDBParser, PPBuilder, Structure
import click


def parse_pdb(pdb: Path,
              structure_id: str,
              suppress_warnings: bool = True) -> Structure.Structure:
    """Parses a PDB file using Bio.PDB and returns the parsed structure.

    Args:
        pdb (Path): Path to the PDB file to be parsed.
        structure_id (str): Identification string for the structure within the PDB file.
        suppress_warnings (bool): Suppress parsing warnings (for efficiency). Defaults to True.

    Returns:
        Structure.Structure: A Bio.PDB Structure object containing the parsed data.

    Raises:
        IOError: If the specified PDB file cannot be found or parsed.
        ValueError: If the structure_id is not found within the PDB file.
    """
    pdb_parser = PDBParser(QUIET=suppress_warnings)

    return pdb_parser.get_structure(structure_id, pdb)

def seq_from_structure(structure: Structure.Structure) -> str:
    """Extracts the amino acid sequence from a Bio.PDB Structure object.

    This function iterates through the chains and peptides within a Bio.PDB Structure
    object and constructs the full amino acid sequence as a string.

    Args:
        structure (Structure.Structure): A Bio.PDB Structure object to extract the sequence from.

    Returns:
        str: The full amino acid sequence as a single string.

    Raises:
        ValueError: If any non-standard amino acid residues are encountered.
    """
    pp_builder = PPBuilder()

    seq_list = []

    for chain in structure.get_chains():

        for pp in pp_builder.build_peptides(chain):

            seq_list.append(str(pp.get_sequence()))

    return ''.join(seq_list)

def close_pipes(pipes: List[Popen]) -> None:
    """Closes the standard output streams of a list of Popen process objects.

    This function facilitates graceful termination of process pipelines by closing
    the standard output streams of the specified processes, allowing them to receive
    SIGPIPE signals if upstream communication terminates.

    Args:
        pipes (List[Popen]): A list of Popen objects representing processes in a pipeline.

    Raises:
        TypeError: If any elements in the pipes list are not Popen objects.
        ValueError: If any of the Popen objects have already terminated.
    """
    for p in pipes:

        p.stdout.close()

def mkdir(context, param, value: Path) -> Path:
    """Creates a directory, ensuring parent directories exist.

    Args:
        context (click.Context): The Click context object for the current command.
            (Not directly used within this function)
        param (click.Parameter): The Click parameter object for the directory option.
            (Not directly used within this function)
        value (Path): The path to the directory to be created.

    Returns:
        Path: The created directory path as a Path object.
    """
    value.mkdir(exist_ok=True, parents=True)

    return value


@click.command()
@click.version_option(version='1.0')
@click.option(
    '--pdb1',
    required=True,
    type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help='PDB 1'
)
@click.option(
    '--pdb2',
    required=True,
    type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help='PDB 2'
)
@click.option(
    '--output_dir',
    required=True,
    type=click.Path(resolve_path=True, path_type=Path),
    callback=mkdir,
    help='Output directory'
)
def cheetah_pdb_format(pdb1: Path, pdb2: Path, output_dir: Path):
    """Cleans and formats two PDB files, extracts amino acid sequences, and creates a combined PDB.

    This function performs the following steps:

    1. Parses the input PDB files using Bio.PDB.
    2. Extracts the amino acid sequences from each parsed structure.
    3. Cleans and formats each PDB file using pdb-tools (https://github.com/haddocking/pdb-tools):
       - Retains only coordinates using pdb_keepcoord.
       - Renumbers residues using pdb_reres.
       - Extracts specified chains using pdb_chain.
       - Tidies the PDB file using pdb_tidy.
    4. Merges the cleaned and formatted PDB files using pdb_merge.
    5. Renumbers residues again and tidies the merged PDB file.
    6. Writes the final combined PDB file to `output_dir/complex_AB.pdb`.
    7. Writes the extracted amino acid sequences to text files in the output directory.

    Args:
        pdb1 (Path): Path to the first PDB file.
        pdb2 (Path): Path to the second PDB file.
        output_dir (Path): Path to the output directory for results.

    Raises:
        IOError: If there are errors reading or writing PDB files.
        subprocess.CalledProcessError: If any `pdb_tools` commands fail.
        ValueError: If unexpected output is encountered from `pdb_tools`.
    """
    chains = ('A', 'B')

    output_pdb = output_dir / 'complex_AB.pdb'

    seq_files = []

    path_in = {chains[0]: pdb1, chains[1]: pdb2}

    path_out = {x: output_dir / f'chain_{x}.pdb' for x in chains}

    for out_chain, pdb in path_in.items():

        struct = parse_pdb(pdb=pdb, structure_id=out_chain)

        # Saved for future if structure chains are needed.
        #struct_chains = [c.get_id() for c in struct.get_chains()]

        seq = seq_from_structure(structure=struct)

        seq_file = output_dir / f'chain_{out_chain}.txt'

        seq_files.append(seq_file)

        with open(seq_file, 'w') as f:

            print(seq, file=f)

        #print(path_out[out_chain], out_chain, sep='\t')

        with open(path_out[out_chain], 'w') as f:

            # Open pipe call and specify commands; final pipe stdout is to file.
            p1 = Popen(['pdb_keepcoord', pdb], stdout=PIPE)

            p2 = Popen(['pdb_reres'], stdin=p1.stdout, stdout=PIPE)

            p3 = Popen(['pdb_chain', f'-{out_chain}'],
                       stdin=p2.stdout, stdout=PIPE)

            p4 = Popen(['pdb_tidy'], stdin=p3.stdout, stdout=f)

            # Allow pipes to receive a SIGPIPE if communication exits.
            close_pipes(pipes=[p1, p2, p3])

            # Execute pipe.
            p4.communicate()

    with open(output_pdb, 'w') as f:

        # Open pipe call and specify commands; final pipe stdout is to file.
        p1 = Popen(['pdb_merge', path_out['A'], path_out['B']], stdout=PIPE)

        p2 = Popen(['pdb_reres'], stdin=p1.stdout, stdout=PIPE)

        p3 = Popen(['pdb_reres'], stdin=p2.stdout, stdout=PIPE)

        p4 = Popen(['pdb_tidy'], stdin=p3.stdout, stdout=f)

        # Allow pipes to receive a SIGPIPE if communication exits.
        close_pipes(pipes=[p1, p2, p3])

        # Execute pipe.
        p4.communicate()


if __name__ == '__main__':

    cheetah_pdb_format()
