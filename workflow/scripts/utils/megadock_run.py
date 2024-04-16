#!/usr/bin/env python3

# TODO:
#   - Include argument to specify the path to the container,
#     and a check to see that the provided container exists.


import subprocess
# from os import getuid
from pathlib import Path

# import docker
from Bio.PDB import PDBParser, Structure


def megadock_run(pdb_a: Path,
                 pdb_b: Path,
                 dock_file: Path,
                 model_number: int, output_dir: Path) -> Structure:
    """
    Runs a Megadock docking simulation using a Singularity container and generates a pseudo-docked PDB structure.

    This function executes a Megadock docking simulation within a Singularity container to generate models for
    two input PDB files (pdb_a and pdb_b). It leverages the 'decoygen' executable within the container to create
    docked models. It then constructs a pseudo-docked PDB by combining the beginning of pdb_a with one of the
    generated decoys and returns the parsed structure.

    Args:
        pdb_a (Path): Path to the first PDB file.
        pdb_b (Path): Path to the second PDB file.
        dock_file (Path): Path to the dock parameter file.
        model_number (int): Number of models to generate.
        output_dir (Path): Path to the output directory for results.

    Returns:
        Structure: A Bio.PDB Structure object representing the pseudo-docked PDB.

    Raises:
        subprocess.CalledProcessError: If the Singularity command or shell commands fail.
        IOError: If there are errors parsing the PDB files.
    """
    # Assign container internal paths:
    # c_data_path = Path('data')

    # c_moving = c_data_path / pdb_b.name

    #c_static = c_data_path / pdb_a.name

    # Build container shell command.
    # megadock_command = ['./decoygen',
    #                     'data/lig_tmp.pdb',
    #                     f'{c_moving}',
    #                     'data/dock.out', f'{model_number}']
    # megadock_command = ['touch', '/opt/MEGADOCK/data/test.txt']

    # Connect to docker daemon.
    # client = docker.from_env()

    # Get current user ID.
    # user_id = getuid()

    # Set MEGADOCK container run command:
    # run = client.containers.run(image="akiyamalab/megadock:cpu-4.1.1",
    #                             command=f"{' '.join(megadock_command)}",
    #                             auto_remove=True,
    #                             runtime='nvidia',
    #                             user=f'{user_id}:{user_id}',
    #                             volumes={output_dir: {
    #                                      'bind': '/opt/MEGADOCK/data',
    #                                      'mode': 'rw'}})

    # Run MEGADOCK container and write STDOUT to log:
    # with open(output_dir / 'log', 'w') as f:

    #     #print(f'\n> Decoy docking {"_"*64}', file=f)
    #     print(run.decode(encoding='utf8'), file=f)

    container = '/srv/data1/home/jo0348st/containers/megadock.sif'

    decoygen = '/usr/local/megadock/megadock-4.1.1/decoygen'

    singularity_command = ['singularity',
                           'exec',
                           '--bind',
                           f'/{str(output_dir).split("/")[1]}',
                           f'{container}',
                           decoygen,
                           f'{output_dir}/lig_tmp.pdb',
                           f'{pdb_b}',
                           f'{dock_file}', f'{model_number}']

    subprocess.run(singularity_command, check=True)

    tmp_lig_pdb = output_dir / 'lig_tmp.pdb'

    tmp_dock_pdb = output_dir / 'dock_tmp.pdb'

    # Build shell command to psuedo-dock the decoy to the other protein.
    os_command_1 = ['head', '-n-1', f'{pdb_a}']

    os_command_2 = ['cat', f'{tmp_lig_pdb}']

    with open(tmp_dock_pdb, 'w') as f:

        subprocess.run(os_command_1, stdout=f)

        subprocess.run(os_command_2, stdout=f)

    # Parse pseudo-docked PDB and return structure object:
    parser = PDBParser()

    return parser.get_structure(f'DOCK_{model_number}', tmp_dock_pdb)
