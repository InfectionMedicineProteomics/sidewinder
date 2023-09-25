import subprocess
# from os import getuid
# from pathlib import Path

# import docker
from Bio.PDB import PDBParser

def megadock_run(pdb_a, pdb_b, dock_file, model_number, output_dir):

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
