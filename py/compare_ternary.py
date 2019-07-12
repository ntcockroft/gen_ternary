#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Compares a list of generated ternary complex PDB files to a reference ternary
 complex PDB file
"""

import argparse
import pandas as pd
from Bio import PDB

def main():
    parser = argparse.ArgumentParser(description='Compares a list of pdb \
                                     files to a reference pdb structure')
    parser.add_argument('-p', '--pdb_file_list', action='store', nargs=1,
                        dest='pdb', help='List of pdb files (.txt) - provide \
                        full pathname if not in current directory')
    parser.add_argument('-r', '--reference_pdb', action='store', nargs=1,
                        dest='ref', help='Full path to reference .pdb file')
    parser.add_argument('-n', '--name', action='store', nargs=1,
                        dest='name', default=['./'], help='Output file name \
                        suffix to be used')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input pdb files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output log should be written')
    args = vars(parser.parse_args())


    #Get list of all pdb files after filtering ZDOCK output
    with open(args['pdb'][0], 'r') as pdb_file_list:
        pdb_files = pdb_file_list.read().splitlines()

    #Get a list of the atom objects from each structure
    parser = PDB.PDBParser()
    structure_atoms = []
    for pdb_file in pdb_files:
        structure = parser.get_structure(pdb_file, args['input'][0] + pdb_file)
        ca_atoms = [atom for atom in structure.get_atoms() if
                    atom.get_id() == 'CA']
        structure_atoms.append(ca_atoms)

    ref_struct = parser.get_structure('ref', args['ref'][0])
    ref_atoms = [atom for atom in ref_struct.get_atoms() if
                 atom.get_id() == 'CA']

    good_pdb = []
    good_rms = []
    all_pdb = []
    all_rms = []
    sup = PDB.Superimposer()
    for idx, struct in enumerate(structure_atoms):
        sup.set_atoms(ref_atoms, struct)
        all_pdb.append(pdb_files[idx])
        all_rms.append(sup.rms)
        if sup.rms <= 10:
            good_pdb.append(pdb_files[idx])
            good_rms.append(sup.rms)

    good_struct = pd.DataFrame()
    good_struct['pdb'] = good_pdb
    good_struct['rms'] = good_rms
    good_struct.to_csv(args['output'][0] + '/' + 'good_'+ args['name'][0] +
                       '.csv', index=False)

    all_struct = pd.DataFrame()
    all_struct['pdb'] = all_pdb
    all_struct['rms'] = all_rms
    all_struct.to_csv(args['output'][0] + '/' + 'all_'+ args['name'][0] +
                      '.csv', index=False)

if __name__ == "__main__":
    main()
