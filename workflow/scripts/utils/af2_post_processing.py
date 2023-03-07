#!/usr/bin/env python3
"""Post-processing of AlphaFold2-Multimer output PDB files.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


from pathlib import Path
from subprocess import Popen, PIPE

# import click


# @click.command()
# @click.version_option(version='1.0')
# @click.option(
#     '--input_pdb',
#     '-i',
#     required=True,
#     type=click.Path(exists=True, path_type=Path),
#     help='AlphaFold2-Multimer output PDB file'
# )
# @click.option(
#     '--chain_pairs',
#     '-c',
#     required=True,
#     type=str,
#     help=(
#         f'Chains to be grouped into output PDBs (format: A,B_C,D'
#         f' where comma separates chains and underscore separates'
#         f' output PDBs)'
#     )
# )
def af2_post_processor(input_pdb: Path, chain_pairs: str):
    """Post-processing of AlphaFold2-Multimer output PDB files.
    Appropriate if AlphaFold2-Multimer was used for docking prior
    to analysis with targeted chemical cross-linking.

    Calls the PDB-tools method pdb_selchain with specified chain
    groupings, pipes them to pdb_reatom to ensure atom numbering
    is correct, and then prints them to new PDB files.
    """

    pdb = input_pdb.resolve()

    # Create list by splitting input string chain pairs on underscore.
    chain_pairs = [x for x in chain_pairs.split('_')]  # Automate with click!

    for chain_pair in chain_pairs:

        # Generate chain pair string "XY" from "X,Y".
        chain_string = ''.join(chain_pair).replace(',', '')

        # Build output file name from input file and chain pair string.
        pdb_out = pdb.parent / f'{pdb.stem}_{chain_string}.pdb'

        with open(pdb_out, 'w') as f:

            # Open pipe call and specify first command (pdb_selchain).
            p1 = Popen(['pdb_selchain', f'-{chain_pair}', pdb], stdout=PIPE)

            # Add next command (pdb_reatom) to pipe and specify output file.
            p2 = Popen(['pdb_reatom'], stdin=p1.stdout, stdout=f)

            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.

            p2.communicate()  # Execute pipe.


if __name__ == '__main__':

    af2_post_processor()
