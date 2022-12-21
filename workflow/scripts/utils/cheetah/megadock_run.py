import subprocess
from os import getuid
from pathlib import Path

import docker
from Bio.PDB import PDBParser

def megadock_run(pdb_a, pdb_b, model_number, output_dir):

    # Assign container internal paths:
    c_data_path = Path('data')

    c_moving = c_data_path / pdb_b.name

    #c_static = c_data_path / pdb_a.name

    # Build container shell command.
    # megadock_command = ['./decoygen',
    #                     'data/lig_tmp.pdb',
    #                     f'{c_moving}',
    #                     'data/dock.out', f'{model_number}']
    megadock_command = ['touch', '/opt/MEGADOCK/data/test.txt']

    # Connect to docker daemon.
    client = docker.from_env()

    # Get current user ID.
    user_id = getuid()

    # Set MEGADOCK container run command:
    run = client.containers.run(image="akiyamalab/megadock:cpu-4.1.1",
                                command=f"{' '.join(megadock_command)}",
                                auto_remove=True,
                                runtime='nvidia',
                                user=f'{user_id}:{user_id}',
                                volumes={output_dir: {
                                         'bind': '/opt/MEGADOCK/data',
                                         'mode': 'rw'}})

    # Run MEGADOCK container and write STDOUT to log:
    # with open(output_dir / 'log', 'w') as f:

    #     #print(f'\n> Decoy docking {"_"*64}', file=f)
    #     print(run.decode(encoding='utf8'), file=f)

    print(run.decode(encoding='utf8'))

    tmp_lig_pdb = output_dir / 'lig_tmp.pdb'

    tmp_dock_pdb = output_dir / 'dock_tmp.pdb'

    print('aopsihdoishgosaihgsaigh')

    # Build shell command to psuedo-dock the decoy to the other protein.
    os_command = ['cat', pdb_a, tmp_lig_pdb, '>', tmp_dock_pdb]

    subprocess.run(os_command)  # Run command.

    # Parse pseudo-docked PDB and return structure object:
    parser = PDBParser()

    return parser.get_structure('DOCK', 'dock_tmp.pdb')
