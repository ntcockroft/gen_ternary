#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Calculates the center of mass distance between the two ligands in the ternary
 complex.
Estimates the hydrophobic solvent-exposed surface area buried in the predicted
 ternary complex.
Filters predicted ternary complexes by ligand center of mass distance  and
 buried hydrophobic surface area.
"""

import argparse
import pymol
import numpy as np

def get_dist(point_1, point_2):
    """
    Obtains the distance between two points

    Args:
        point_1: a numpy array of coordinates
        point_2: a numpy array of coordinates

    Returns:
        A numpy.float64 corresponding to the distance between the two given
        points
    """
    squared_dist = np.sum((point_1-point_2)**2, axis=0)
    dist = np.sqrt(squared_dist)

    return dist


def get_hydrophobic_sasa(pdb_file):
    """
    Approximates hydrophobic solvent-accessible surface area by calculating
    the solvent-accessible surface area of hydrophobic residues.

    Args:
        pdb_file: Path to pdb file for calculation

    Returns:
        Approximate hydrophobic solvent-accessible surface area of the pdb file
    """
    # use solvent-accessible surface with high sampling density
    pymol.cmd.set('dot_solvent', 1)
    pymol.cmd.set('dot_density', 3)

    pdb_file_name = pdb_file.split('/')[-1].split('.')[0]
    pymol.cmd.load(pdb_file)
    pymol.cmd.create(pdb_file_name, 'all')
    obj_list = pymol.cmd.get_object_list('all')
    if len(obj_list) > 1:
        for obj in obj_list[0:-1]:
            pymol.cmd.delete(obj)
    pymol.cmd.remove('hetatm')
    pymol.cmd.select('hydrophobes', '(resn ala+gly+val+ile+leu+phe+met)')
    hydrophobic_sasa = pymol.cmd.get_area('hydrophobes')
    pymol.cmd.delete(pdb_file_name)
    return hydrophobic_sasa


def main():
    parser = argparse.ArgumentParser(description='Filters protein-protein pdb \
                                     structures that were output by Z-dock')
    parser.add_argument('-p', '--pdb_file', action='store', nargs=1,
                        dest='pdb', help='Z-dock result file (.pdb)')
    parser.add_argument('-a', '--pdb_a', action='store', nargs=1,
                        dest='pdb_a', help='Reference pdb of a protein in the \
                        complex (.pdb) - use full pathname')
    parser.add_argument('-b', '--pdb_b', action='store', nargs=1,
                        dest='pdb_b', help='Reference pdb of a protein in the \
                        complex (.pdb) - use full pathname')
    parser.add_argument('-l', '--ligand_names', action='store', nargs=2,
                        dest='ligand', help='Names of ligands present in each \
                        docked protein')
    parser.add_argument('-d', '--distance', action='store', nargs=1,
                        type=float, dest='dist', help='Max distance between \
                        the center of mass of each ligand - should be \
                        approximately the length of the fully elongated \
                        PROTAC')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input data files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output files should be written')
    args = vars(parser.parse_args())

    #Calculate the hydrophobic solvent-accessible surface area buried by the
    #predicted complex
    complex_sasa = get_hydrophobic_sasa(args['input'][0] + '/'
                                        + args['pdb'][0])
    pdb_a_sasa = get_hydrophobic_sasa(args['pdb_a'][0])
    pdb_b_sasa = get_hydrophobic_sasa(args['pdb_b'][0])

    buried_sasa = (pdb_a_sasa + pdb_b_sasa) - complex_sasa

    #Calculate center of mass distance between ligands
    pymol.cmd.load(args['input'][0] + '/' + args['pdb'][0])
    pymol.cmd.select('lig_1', 'resn ' + args['ligand'][0])
    pymol.cmd.select('lig_2', 'resn ' + args['ligand'][1])
    lig_com_1 = pymol.cmd.centerofmass('lig_1')
    lig_com_2 = pymol.cmd.centerofmass('lig_2')
    lig_com_1 = np.array(lig_com_1)
    lig_com_2 = np.array(lig_com_2)

    dist = get_dist(lig_com_1, lig_com_2)


    if dist <= args['dist'][0] and buried_sasa >= 100:
        pymol.cmd.save(args['output'][0] + '/' + args['pdb'][0])
        print(args['pdb'][0] + ' ' + str(dist))


if __name__ == "__main__":
    main()
