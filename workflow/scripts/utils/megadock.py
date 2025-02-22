#!/usr/bin/env python3
"""Python docker SDK code calling MEGADOCK 4.

[2023-01-04] Currently not functioning: discontinuous memory issue.
"""

__author__ = 'Joel Ströbaek'
__email__ = 'joel.strobaek@gmail.com'


from os import getuid
from pathlib import Path

import docker
import click

@click.command()
@click.version_option(version='1.0')
@click.option('--receptor',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--ligand',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--output_dir',
              required=True,
              type=click.Path(exists=True, path_type=Path), help='')
@click.option('--n_out',
              type=int,
              default=10000,
              help='Number of output models')
def run_megadock(receptor: Path,
                 ligand: Path, output_dir: Path, n_out: int) -> Path:
    """Dock two proteins using a containerized MEGADOCK 4.

    ...

    Parameters
    ----------

    """

    # Assign container internal paths:
    c_data_path = Path('data')

    c_receptor = c_data_path / receptor.name

    c_ligand = c_data_path / ligand.name

    c_docked_out = c_data_path / f'{receptor.stem}-{ligand.stem}.out'

    # Build container shell command.
    megadock_command = ['megadock-gpu',
                        f'-R {c_receptor}',
                        f'-L {c_ligand}',
                        f'-N {n_out}',
                        f'-o {c_docked_out}']

    # Connect to docker daemon.
    client = docker.from_env()

    # Get current user ID.
    user_id = getuid()

    # Set MEGADOCK container run command:
    run = client.containers.run(image="akiyamalab/megadock:gpu-4.1.4",
                                command=' '.join(megadock_command),
                                auto_remove=True,
                                runtime='nvidia',
                                user=user_id,
                                volumes={output_dir: {'bind':
                                                      '/opt/MEGADOCK/data',
                                                      'mode':
                                                      'rw'}})

    # Run MEGADOCK container and write STDOUT to log:
    with open(output_dir / 'log', 'w') as f:

        print(run.decode(encoding='utf8'), file=f)

    #return output_dir / c_docked_out.name

if __name__ == '__main__':

    run_megadock()
