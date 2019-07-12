#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Re-orders the chains in a PDB file based on a user specified ordering. ZRANK
 compares the first two chains in a PDB file and if a ternary complex has more
 than 2 chains present, the order may need to be modified in the PDB file to
 allow for correct scoring.
"""

import argparse
from Bio import PDB

def main():
    parser = argparse.ArgumentParser(description='Re-orders chains in the PDB \
                                     file based on the specified ordering')
    parser.add_argument('-p', '--pdb', action='store', nargs=1,
                        dest='pdb', help='PDB file with chains to be ordered')
    parser.add_argument('-r', '--rank_order', action='store', nargs='*',
                        dest='order', type=int, help='Space separated rank \
                        order for chains. Top rank starts with 0 \
                        (instead of 1)')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input pdb files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output log should be written')
    args = vars(parser.parse_args())

    parser = PDB.PDBParser()
    structure = parser.get_structure('pdb',
                                     args['input'][0] + '/' + args['pdb'][0])

    chain_list = []
    for model in structure:
        for chain in model:
            chain_list.append(chain)

    for chain in chain_list:
        model.detach_child(chain.id)

    new_order = args['order']
    chain_list = [chain_list[i] for i in new_order]

    for idx, chain in enumerate(chain_list):
        model.insert(idx, chain)

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(args['output'][0] + '/' + args['pdb'][0])

if __name__ == "__main__":
    main()
