#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ntcockroft

Generates conformers of a PROTAC molecule and identifies conformers that could
 reasonably span the binding sites of the ligase and target protein structures.
Many conformations can be generated, but only a few should be needed to get
 close matches due to the use of positional restraints.
"""

import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
import pybel

def read_pdb(pdb_file):
    """
    Reads in a .pdb file and returns and RDKit Mol object. The RDKit
     Chem.MolFromPDBFile function loses bond order information. The Pybel
     readfile function reatains bond order information. Ultimately, the RDKit
     Mol object is produced by leveraging these two packages.

    Args:
        pdb_file: A small molecule .pdb file

    Returns:
        pdb_mol_h_bo: A RDKit Mol object for the molecule with hydrogen atoms
                      and bond orders
    """

    for molecule in pybel.readfile('pdb', pdb_file):
        ref_pdb_smi = molecule.write('smiles').split()[0]
        ref_pdb_mol = Chem.MolFromSmiles(ref_pdb_smi)
        ref_pdb_mol_h = Chem.AddHs(ref_pdb_mol)

    pdb_mol_h = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    pdb_mol_h_bo = AllChem.AssignBondOrdersFromTemplate(ref_pdb_mol_h,
                                                        pdb_mol_h)

    return pdb_mol_h_bo


def mcs_atom_map(prb_mol, ref_mol):
    """
    Obtains and atom map for the maximum common substructure of two molecules

    Args:
        prb_mol: RDKit Mol object of the probe molecule
        ref_mol: RDKit Mol object of the reference molecule

    Returns:
        atom_maps: A list that maps the atoms from the probe molecule to the
                  atoms of the reference molecule in their maximum common
                  substructure
        pattern: The RDKit Mol object of the maximum common substructure
    """

    mols = [prb_mol, ref_mol]
    mcs_result = rdFMCS.FindMCS(mols)
    pattern = Chem.MolFromSmarts(mcs_result.smartsString)

    prb_atoms = prb_mol.GetSubstructMatches(pattern)
    ref_atoms = ref_mol.GetSubstructMatch(pattern)

    atom_map = list(zip(prb_atoms[0], ref_atoms))

    return atom_map, pattern


def restr_coords_mcs(prb_mol, ref_mol, refCid=-1):
    """
    Maps the atoms of the probe molecule to the coordinates of matching atoms
     in the reference molecule. The coordinates are intended for use as
     positional restraints in conformer geneartion of the probe molecule.

    Args:
        prb_mol: RDKit Mol object of the probe molecule
        ref_mol: RDKit Mol object of the reference molecule

    Returns:
        atom_maps: A list that maps the atoms from the probe molecule to the
                  atoms of the reference molecule in their maximum common
                  substructure
        pattern: The RDKit Mol object of the maximum common substructure
    """

    atom_map, pattern = mcs_atom_map(prb_mol, ref_mol)
    ref_conf = ref_mol.GetConformer(refCid)

    restr_coords = {}
    for pair in atom_map:
        restr_coords[pair[0]] = ref_conf.GetAtomPosition(pair[1])

    return restr_coords


def rmsd(prb_mol, ref_mol, prbCid=-1, refCid=-1):
    """
    Obtains the root mean squared distance between the maximum common
     substructure of two molecules

    Args:
        prb_mol: RDKit Mol object of the probe molecule conformers
        ref_mol: RDKit Mol object of the reference molecule conformers
        prbCid: The conformer id of the desired probe molecule conformer
        refCid: The conformer id of the desired probe molecule conformer

    Returns:
        The distance as a float between the two molecule's maximum common
         substructure
    """

    atom_map, pattern = mcs_atom_map(prb_mol, ref_mol)
    prb_conf = prb_mol.GetConformer(prbCid)
    ref_conf = ref_mol.GetConformer(refCid)

    prb_coords = []
    ref_coords = []
    for pair in atom_map:
        prb_pt = prb_conf.GetAtomPosition(pair[0])
        ref_pt = ref_conf.GetAtomPosition(pair[1])
        prb_coords.append([prb_pt.x, prb_pt.y, prb_pt.z])
        ref_coords.append([ref_pt.x, ref_pt.y, ref_pt.z])

    dim = len(prb_coords[0])
    n_points = len(prb_coords)
    result = 0.0
    for v, w in zip(prb_coords, ref_coords):
        result += sum([(v[i] - w[i])**2.0 for i in range(dim)])
    dist = np.sqrt(result/n_points)

    return dist


def main():
    parser = argparse.ArgumentParser(description='Generates PROTAC conformers \
                                     and filters based on distance to complex \
                                     ligands.')
    parser.add_argument('-p', '--protac_smi', action='store', nargs=1,
                        dest='protac', help='SMILES string of PROTAC ligand')
    parser.add_argument('-t', '--target_ligand', action='store', nargs=1,
                        dest='target', help='Target protein ligand (.pdb)')
    parser.add_argument('-l', '--ligase_ligand', action='store', nargs=1,
                        dest='ligase', help='Ligase protein ligand')
    parser.add_argument('-n', '--num_conf', action='store', nargs=1, type=int,
                        dest='num', default=[10], help='Number of conformers \
                        to generate')
    parser.add_argument('-r', '--rmsd', action='store', nargs=1, type=float,
                        dest='rmsd', help='RMSD threshold for PROTAC \
                        conformer')
    parser.add_argument('-c', '--num_cores', action='store', nargs=1, type=int,
                        dest='cores', default=[0], help='Number of cores to \
                        use for conformer generation')
    parser.add_argument('-i', '--input_directory', action='store', nargs=1,
                        dest='input', default=['./'], help='Directory where \
                        input data files are stored')
    parser.add_argument('-o', '--output_directory', action='store', nargs=1,
                        dest='output', default=['./'], help='Directory where \
                        output files should be written')
    args = vars(parser.parse_args())

    #Get SMILES of PROTAC structure and add hydrogens
    protac = Chem.MolFromSmiles(args['protac'][0])
    protac = Chem.AddHs(protac)

    #Get ligands from protein-protein docking result - ligands have to be
    #extracted and saved as individaul .pdb files
    target_lig = read_pdb(args['input'][0] + '/' + args['target'][0])
    ligase_lig = read_pdb(args['input'][0] + '/' +  args['ligase'][0])

    #Obtain restraints for the protac molecule
    protac_ligase_coords = restr_coords_mcs(protac, ligase_lig)
    protac_target_coords = restr_coords_mcs(protac, target_lig)
    protac_coords = {**protac_ligase_coords, **protac_target_coords}


    #Generate conformations with positional restraints
    cids = AllChem.EmbedMultipleConfs(protac, numConfs=args['num'][0],
                                      numThreads=args['cores'][0],
                                      coordMap=protac_coords)

    """
    Cycle through all generated protac conformers, aligns the ligase binding
    portion, and then compares the distance between the target protein ligand
    from the protein-protein docking result and the analogous portion of the
    PROTAC.
    """
    atom_map, mcs = mcs_atom_map(protac, ligase_lig)
    for i in cids:
        #Alignment without reflection
        Chem.rdMolAlign.AlignMol(prbMol=protac, refMol=ligase_lig,
                                 prbCid=i, atomMap=atom_map)
        rms_1 = rmsd(protac, target_lig, prbCid=i)

        #Alignmnet with reflection
        Chem.rdMolAlign.AlignMol(prbMol=protac, refMol=ligase_lig,
                                 prbCid=i, atomMap=atom_map, reflect=True)
        rms_2 = rmsd(protac, target_lig, prbCid=i)

        #Choos alignment with best RMSD
        if rms_1 > rms_2:
            Chem.rdMolAlign.AlignMol(prbMol=protac, refMol=ligase_lig,
                                     prbCid=i, atomMap=atom_map)
            rms = rmsd(protac, target_lig, prbCid=i)
        else:
             Chem.rdMolAlign.AlignMol(prbMol=protac, refMol=ligase_lig,
                                      prbCid=i, atomMap=atom_map, reflect=True)
             rms = rmsd(protac, target_lig, prbCid=i)

        if i == 0:
            print('complex.' + args['output'][0].split('/')[-1] + '.pdb')
        fname = args['output'][0] + '/' + 'protac_conf_' + str(i) + '.pdb'
        Chem.MolToPDBFile(protac, filename=fname, confId=i)
        print('protac_conf_' + str(i) + '.pdb', rms)


if __name__ == "__main__":
    main()
