#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Performs a local optimization of a small molecule ligand using MMFF94 in
 pybel (openbabel)
"""

import argparse
import pybel

def local_opt_ligand(ligand_pdb, output_dir):
    """
    Locally optimizes the coordinates of the provided small molecule. Uses
    the MMFF94 forcefield by default.

    Args:
        ligand_pdb: PDB file of the small molecule to optimize
        output_dir: Directory to write optimized small molecules to
    Returns:
        Nothing. Writes and SDF and PDB file of the optimize
        ligand to the specified output directory.
    """
    lig_name = ligand_pdb.split('.pdb')[-2].split('/')[-1]
    for molecule in pybel.readfile('pdb', ligand_pdb):
        molecule.localopt(steps=1000)
        molecule.write('sdf', output_dir + '/' + lig_name +
                       '_opt.sdf', overwrite=True)
        molecule.write('pdb', output_dir + '/' + lig_name +
                       '_opt.pdb', overwrite=True)


def main():
    parser = argparse.ArgumentParser(description='Runs a local optimization \
                                     on a provided small molecule pdb file')
    parser.add_argument('-l', '--ligand', action='store', nargs=1,
                        dest='ligand', help='The ligand .pdb file to generate \
                        parameters for')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input pdb files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output log should be written')
    args = vars(parser.parse_args())

    lig_pdb = args['input'][0] + '/' + args['ligand'][0]
    local_opt_ligand(lig_pdb, args['output'][0])


if __name__ == "__main__":
    main()
