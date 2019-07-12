#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Performs a local minimzation on predicted ternary complex structures using
 OpenMM
Prior to running this minimization, ligand parameters need to be generated
 using the param_lig.py script.
"""

import argparse
import parmed
from parmed.openmm import load_topology
from simtk.openmm import LangevinIntegrator
from simtk.openmm.app import PDBFile, ForceField, Simulation, HBonds, NoCutoff
from simtk.unit import angstrom, kelvin, picoseconds, picosecond

def main():
    parser = argparse.ArgumentParser(description='Runs local minimization on \
                                     the protein-ligand complex using OpenMM')
    parser.add_argument('-l', '--ligand', action='store', nargs=1,
                        dest='ligand', help='The ligand .pdb file to generate \
                        parameters for')
    parser.add_argument('-x', '--xml_ligand', action='store', nargs=1,
                        dest='xml', help='The xml file containing ligand \
                        parameters - use full path if not in directory where \
                        this script is being executed from')
    parser.add_argument('-p', '--protein', action='store', nargs=1,
                        dest='protein', help='The protein complex .pdb file')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input pdb files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output log should be written')
    args = vars(parser.parse_args())

    #Load in ligand file and parameters
    ligand_pdbfile = PDBFile(args['input'][0] + '/' + args['ligand'][0])
    ligand_system = parmed.load_file(args['xml'][0])
    ligand_structure = load_topology(ligand_pdbfile.topology, ligand_system,
                                     xyz=ligand_pdbfile.positions)

    #Load in protein complex file and  force field
    receptor_pdbfile = PDBFile(args['input'][0] + '/' + args['protein'][0])
    receptor_pdbfile_name = args['protein'][0].split('.pdb')[-2]
    omm_forcefield = ForceField('amber14-all.xml')
    receptor_system = omm_forcefield.createSystem(receptor_pdbfile.topology)
    receptor_structure = load_topology(receptor_pdbfile.topology,
                                       receptor_system,
                                       xyz=receptor_pdbfile.positions)

    #Generate ligand-protein complex
    complex_structure = receptor_structure + ligand_structure
    complex_system = (complex_structure
                      .createSystem(nonbondedMethod=NoCutoff,
                                    nonbondedCutoff=9.0*angstrom,
                                    constraints=HBonds,
                                    removeCMMotion=False))

    #Set up simulation parameters
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond,
                                    0.002*picoseconds)
    simulation = Simulation(complex_structure.topology, complex_system,
                            integrator)
    simulation.context.setPositions(complex_structure.positions)

    #Run local minimization
    simulation.minimizeEnergy()
    positions = simulation.context.getState(getPositions=True).getPositions()

    #Save minimized complex structure
    PDBFile.writeFile(simulation.topology, positions,
                      open(args['output'][0] + '/' + receptor_pdbfile_name +
                           '_min.pdb', 'w'))


if __name__ == "__main__":
    main()
