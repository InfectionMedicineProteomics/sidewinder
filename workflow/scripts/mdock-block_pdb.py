#!/usr/bin/env python3

import copy
import shutil
import warnings
from dataclasses import dataclass
from pathlib import Path
from string import ascii_uppercase
from typing import List

import click
from abnumber import Chain
from Bio import BiopythonParserWarning, SeqIO
from typing_extensions import Literal

warnings.simplefilter('ignore', BiopythonParserWarning)


@dataclass
class block_pdb:
    """Class for blocking the input receptor for MegaDock"""
    receptor: Path

    output_dir: Path = Path.cwd() / 'megadock_out'

    chain_ids = [*ascii_uppercase]

    def __post_init__(self):

        # def get_chains(name: str, pdb: Path):

        #     parser = PDBParser()

        #     structure = parser.get_structure(name, pdb)

        #     chains = [x.get_id() for x in structure.get_chains()]

        #     return chains

        self.output_dir.mkdir(exist_ok=True)

        self.receptor_blocked = None

        # self.receptor_chains = get_chains(self.receptor.stem, self.receptor)

        # self.ligand_chains = get_chains(self.ligand.stem, self.ligand)

    def get_resarg(self, res_index: List[int]) -> str:

        n_seq = []

        res_arg_list = []

        for i in range(len(res_index)):

            n = res_index[i]

            try:

                n_next = res_index[i+1]

            except:

                if len(n_seq) > 0:

                    res_arg_list.append(f'{n_seq.pop()}-{n}')

                else:

                    res_arg_list.append(f'{n}')

                break

            if n+1 == n_next and len(n_seq) == 0:

                n_seq.append(n)

            elif n+1 != n_next and len(n_seq) == 1:

                res_arg_list.append(f'{n_seq.pop()}-{n}')

            elif len(n_seq) == 0:

                res_arg_list.append(str(n))

        return ','.join(res_arg_list)

    def block_receptor(self, chain: str, res_block: str):

        unlink_prev = False

        if self.receptor_blocked is not None:

            pdb = self.receptor_blocked

            output_file = self.output_dir / f'{pdb.stem}_{chain}.pdb'

            unlink_prev = True

        else:

            pdb = self.receptor

            output_file = self.output_dir / f'{pdb.stem}_block_{chain}.pdb'

        res_list = []

        for r in res_block.split(','):

            if '-' in r:

                r_min, r_max = r.split('-')

                i_range = [str(i) for i in range(int(r_min), int(r_max) + 1)]

                res_list.extend(i_range)

            else:

                res_list.append(r)

        with open(pdb, 'r') as f_in, open(output_file, 'w') as f_out:

            for l in f_in.readlines():

                l = l.strip()

                if l[0:4] not in ['ATOM', 'TER '] and l[0:6] != 'HETATM':

                    f_out.write(f'{l}\n')

                    continue

                if l[21] != chain:

                    f_out.write(f'{l}\n')

                    continue

                ll = l

                if ll[22:26].strip() not in res_list:

                    f_out.write(f'{l}\n')

                    continue

                f_out.write(f'{l[0:16]} BLK{l[20:]}\n')

        if unlink_prev:

            self.receptor_blocked.unlink()

        self.receptor_blocked = Path(output_file)


@dataclass
class fv_indexing:
    """..."""

    fv: Path

    scheme: Literal['chothia', 'imgt', 'kabat'] = 'chothia'

    def __post_init__(self):

        fv_dict = {'CDR1': [], 'CDR2': [], 'CDR3': [],
                   'FR1': [], 'FR2': [], 'FR3': [], 'FR4': []}

        self.seq_records = {}

        self.chains = {}

        chain_types = {}

        for rec in SeqIO.parse(self.fv, 'pdb-atom'):

            chain = rec.annotations['chain']

            seq = str(rec.seq)

            self.seq_records[chain] = seq

            self.chains[chain] = Chain(seq, scheme=self.scheme, name=chain)

            chain_types[chain] = self.chains[chain].chain_type

        chain_annots = {chain: {'chain_type': chain_types[chain],
                                'index': copy.deepcopy(fv_dict)}
                        for chain in [*self.chains.keys()]}

        for chain in self.chains.values():

            idx = 1

            for pos, res in chain:

                region = pos.get_region()

                set_index = chain_annots[chain.name]['index'][region]

                set_index.append(idx)

                idx += 1

        self.chain_annotations = chain_annots


@click.command()
@click.version_option(version='1.0a')
@click.option('--multi_chain_pdb',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='Path to Sidewinder-MS input antibody Fv .pdb')
@click.option('--single_chain_pdb',
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='Path to Sidewinder-MS edited single chain antibody Fv .pdb')
@click.option('--output_dir',
              '-o',
              required=True,
              type=click.Path(path_type=Path),
              help='Path to output directory')
def block_fv_pdb(multi_chain_pdb, single_chain_pdb, output_dir):
    """..."""

    multi_chain_pdb = multi_chain_pdb.resolve()

    single_chain_pdb = single_chain_pdb.resolve()

    output_dir = output_dir.resolve()

    fv_index = fv_indexing(multi_chain_pdb, 'chothia')

    blocker = block_pdb(single_chain_pdb, output_dir)

    index_convert = {}

    last_index = 1

    for chain, seq in fv_index.seq_records.items():

        index_convert[chain] = {}

        for idx, res in enumerate(seq):

            index_convert[chain][idx+1] = last_index

            last_index += 1

    block_list = []

    for chain, d in fv_index.chain_annotations.items():

        block_index = []

        for region, index in d['index'].items():

            # Blocks the FR of the Fv region.
            if 'CDR' not in region[:-1]:

                block = [index_convert[chain][idx] for idx in index]

                block_index.extend(block)

        # Also block the Fc region of current chain.
        block = [index_convert[chain][idx + 1]
                 for idx in range(len(fv_index.chains[chain].seq),
                                  len(fv_index.seq_records[chain]))]

        block_index.extend(block)

        block_list.append(blocker.get_resarg(block_index))

    block_string = ','.join(block_list)

    blocker.block_receptor('A', block_string)

    tmp_out = blocker.receptor_blocked

    final_out = str(output_dir / 'chain_A_blocked.pdb')

    shutil.copy(str(tmp_out), final_out)

    tmp_out.unlink()


if __name__ == "__main__":

    block_fv_pdb()
