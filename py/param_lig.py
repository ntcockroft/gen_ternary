#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Parameterizes a small molecule ligand using the SMIRNOFF99Frosst forcefield
 from openforcefield for use in a OpenMM minimizaton
Prior to parameterizing the small molecule ligand, the geometry should be
 fixed using the minimize_lig.py script. The restrained RDKit conformer
 generation can produce molecules that have incorrect bond angles/distances due
 to the restraints, which need to be resolved before paramterization.
"""

import argparse
import time
from simtk.openmm import XmlSerializer
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField

def main():
    parser = argparse.ArgumentParser(description='Parameterizes a small \
                                     molecule ligand for use with OpenMM \
                                     using OpenFF')
    parser.add_argument('-l', '--ligand', action='store', nargs=1,
                        dest='ligand', help='The ligand .sdf file to generate \
                        parameters for')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input pdb files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output log should be written')
    args = vars(parser.parse_args())

    #Load SDF file from minimize_lig.py
    lig_sdf = args['input'][0] + '/' + args['ligand'][0]
    lig_name = lig_sdf.split('.sdf')[-2]

    lig_off_molecule = Molecule(args['output'][0] + '/' + lig_sdf)
    force_field = ForceField('test_forcefields/smirnoff99Frosst.offxml')
    start = time.time()
    ligand_system = force_field.create_openmm_system(lig_off_molecule
                                                     .to_topology())
    end = time.time()
    print(end - start)

    with open(lig_name + '.xml', 'w') as f:
        f.write(XmlSerializer.serialize(ligand_system))


if __name__ == "__main__":
    main()
