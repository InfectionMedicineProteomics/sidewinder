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
from typing import Tuple
import warnings

from Bio.PDB import PDBParser, PPBuilder, Structure
import click


def parse_pdb(pdb: Path, structure_id: str) -> Structure.Structure:
    """Get structure from specified PDB file.
    """

    pdb_parser = PDBParser(QUIET=True)

    return pdb_parser.get_structure(structure_id, pdb)

def seq_from_structure(structure: Structure.Structure) -> str:
    """Get sequence from structure object.
    """

    pp_builder = PPBuilder()

    seq_list = []

    for chain in structure.get_chains():

        for pp in pp_builder.build_peptides(chain):

            seq_list.append(str(pp.get_sequence()))

    return ''.join(seq_list)

def close_pipes(pipes: list) -> None:
    """Close each listed pipe.
    """

    for p in pipes:

        p.stdout.close()


@click.command()
@click.version_option(version='1.0')
@click.option(
    '--pdb1',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='PDB 1'
)
@click.option(
    '--pdb2',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='PDB 2'
)
@click.option(
    '--output_dir',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Output directory'
)
def cheetah_pdb_format(pdb1: Path,
    pdb2: Path, output_dir: Path) -> Tuple[Path, Path, Path]:
    """
    """

    chains = ('A', 'B')

    output_pdb = output_dir / 'complex_AB.pdb'

    seq_files = []

    path_in = {chains[0]: pdb1, chains[1]: pdb2}

    path_out = {x: output_dir / f'model_{x}.pdb' for x in chains}

    for out_chain, pdb in path_in.items():

        struct = parse_pdb(pdb=pdb, structure_id=out_chain)

        # Saved for future if structure chains are needed.
        #struct_chains = [c.get_id() for c in struct.get_chains()]

        seq = seq_from_structure(structure=struct)

        seq_file = output_dir / f'seq_{out_chain}.txt'

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

    #return output_pdb, seq_files[0], seq_files[1]

if __name__ == '__main__':

    cheetah_pdb_format()
