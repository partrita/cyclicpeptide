#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File Name: structure2sequence.py
Author: Dingfeng Wu
Creator: Dingfeng Wu
Date Created: 2022-11-15
Last Modified: 2023-12-26
Version: 1.0.1
License: MIT License
Description: Structure-to-Sequence (Struc2seq) is a computing process based on RDkit and the characteristics of cyclic peptide sequences, which can convert cyclic peptide SMILES into sequence information.

Copyright Information: Copyright (c) 2023 dfwlab (https://dfwlab.github.io/)

The code in this script can be used under the MIT License.
"""

from io import BytesIO
from itertools import product
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import SVG, display
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import pkg_resources
from .setting import *
from rdkit.Chem.rdchem import RWMol, AtomPDBResidueInfo
def mol2seq_for_essentialAA(m):
    """
    **Using the RDKit function "Chem.MolToSequence" conversion**

    * Can convert peptide sequences with essential amino acids;
    * Sequences containing modified or non-essential amino acids cannot be accurately identified;
    * Suitable for linear peptides, unable to recognize the cyclization position of cyclic peptides;

    Example::

        smiles = 'CCC(C)[C@H](N)C1=NCC(C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CCCCNC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CC(=O)O)NC(=O)[C@H](Cc3cnc[nH]3)NC(=O)[C@@H](Cc3ccccc3)NC(=O)[C@H](C(C)CC)NC(=O)[C@@H](CCCN)NC2=O)C(C)CC)S1'
        m = Chem.MolFromSmiles(smiles)
        m_renum, seq = struc2seq.mol2seq_for_essentialAA(m)
        print('SEQ:', seq)

    SEQ: IRHTCCVGVCFLMACICIEQFDPCEM
    """
    aa_smiles = {}
    file_path1 = aa_smiles_path  # Replace with your file path
    print(file_path1)
    with open(file_path1, 'r') as file:
        for line in file:
            if line.strip():  # Skip blank lines
                aa, smile = line.strip().split(':')
                aa_smiles[aa.strip()] = smile.strip()
    file_path2 = aas_path  # Replace with your file path
    with open(file_path2, 'r') as file:
        aas = [line.strip() for line in file if line.strip()]
        # order important because gly is substructure of other aas
    # detect the atoms of the backbone and assign them with info
    CAatoms = m.GetSubstructMatches(Chem.MolFromSmarts("[C:0](=[O:1])[C:2][N:3]"))
    # print(CAatoms)
    for atoms in CAatoms:
        a = m.GetAtomWithIdx(atoms[2])
        info = Chem.AtomPDBResidueInfo()
        info.SetName(" CA ")  # spaces are important
        a.SetMonomerInfo(info)
    # detect the presence of residues and set residue name for CA atoms only
    for curr_aa in aas:
        matches = m.GetSubstructMatches(Chem.MolFromSmiles(aa_smiles[curr_aa]))
        for atoms in matches:
            for atom in atoms:
                a = m.GetAtomWithIdx(atom)
                info = Chem.AtomPDBResidueInfo()
                if a.GetMonomerInfo() is not None:
                    if a.GetMonomerInfo().GetName() == " CA ":
                        info.SetName(" CA ")
                        info.SetResidueName(curr_aa)
                        a.SetMonomerInfo(info)
    # renumber the backbone atoms so the sequence order is correct:
    # N-terminal acetylation and C-terminal amidation: https://zhuanlan.zhihu.com/p/540769025
    # Peptide chains are read from the N-terminus to the C-terminus; the C-terminus usually keeps a COOH (HO-C(=O)), followed by one C and then peptide bonds NC(=O) all the way to the N-terminus
    # Both forward and reverse recognition work; reversed it is OC(=O)CNC(=O)CNC(=O)CN..., or NC(=O)CNC(=O)CNC(=O)CN... if the C-terminus is amidated
    # The N-terminus always keeps N: -NH2 without acetylation, -NC(=O)CH3 with acetylation;
    bbsmiles = "C(=O)CN" * len(m.GetSubstructMatches(Chem.MolFromSmiles(aa_smiles["GLY"])))  # generate backbone SMILES
    #backbone = m.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))[0]
    matches = m.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))
    backbone = ''
    if matches:  # Check whether there are matches
        backbone = matches[0]
        # Other operations
    else:
        print("No matches found.")
    #if not len(backbone):
        #return "Not AA"
    id_list = list(backbone)
    id_list.reverse()
    # print(id_list)
    for idx in [a.GetIdx() for a in m.GetAtoms()]:
        if idx not in id_list:
            id_list.append(idx)
    # print(id_list)
    m_renum = Chem.RenumberAtoms(m, newOrder=id_list)
    return m_renum, Chem.MolToSequence(m_renum)

def blend_colors(c1, c2, alpha):
    """
    Blend two colors based on a given alpha value.

    :param c1: The first color represented as a tuple of three values (red, green, blue).
    :type c1: tuple[int, int, int]
    :param c2: The second color represented as a tuple of three values (red, green, blue).
    :type c2: tuple[int, int, int]
    :param alpha: The blending factor, a value between 0 and 1 that determines the proportion of each color in the blend.
    :type alpha: float
    :return: A tuple representing the blended color with values for red, green, and blue.
    :rtype: tuple[int, int, int]

    This function takes two colors (`c1` and `c2`) represented as tuples of red, green, and blue values, along with a blending factor (`alpha`). It then calculates the blended color by applying a weighted average formula for each color channel (red, green, and blue). The resulting blended color is returned as a tuple of the new red, green, and blue values.
    """

    # Unpack the colors
    r1, g1, b1 = c1
    r2, g2, b2 = c2
    # Compute the blended color
    r = r1 * alpha + r2 * alpha * (1 - alpha)
    g = g1 * alpha + g2 * alpha * (1 - alpha)
    b = b1 * alpha + b2 * alpha * (1 - alpha)
    return r, g, b

def detect_backbone(mol):
    """
    Detect the backbone of a given molecule.

    :param mol: The molecule object for which the backbone is to be detected.
    :type mol: Chem.Mol
    :return: A tuple containing two elements:
        - The first element is the backbone of the molecule represented as a tuple of atom indices if matches are found; otherwise, an empty string.
        - The second element is a list of the backbone atom indices in reverse order if matches are found; otherwise, None.
    :rtype: tuple[str, list[int] or None]

    This function aims to identify the backbone of a molecule. 
    """
    # Identify the backbone and number backbone atoms 0..N; side-chain atom indices are larger than N
    # renumber the backbone atoms so the sequence order is correct:
    # N-terminal acetylation and C-terminal amidation: https://zhuanlan.zhihu.com/p/540769025
    # Peptide chains are read from the N-terminus to the C-terminus; the C-terminus usually keeps a COOH (HO-C(=O)), followed by one C and then peptide bonds NC(=O) all the way to the N-terminus
    # Both forward and reverse recognition work; reversed it is OC(=O)CNC(=O)CNC(=O)CN..., or NC(=O)CNC(=O)CNC(=O)CN... if the C-terminus is amidated
    # The N-terminus always keeps N: -NH2 without acetylation, -NC(=O)CH3 with acetylation;
    bbsmiles = "C(=O)CN" * len(
        mol.GetSubstructMatches(Chem.MolFromSmiles('NCC=O')))  # Count residues with the minimal GLY unit and generate the backbone SMILES
    # backbone = mol.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))[0] raises "tuple index out of range"
    matches = mol.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))
    backbone = ''
    if matches:  # Check whether there are matches
        backbone = matches[0]
        # Other operations
    else:
        print("No matches found.")
        return backbone, None # Exit if no peptide bonds are found or the counts differ, e.g., two linear chains cyclized via two links must be searched one by one;
    backbone_idx = list(backbone)
    backbone_idx.reverse()
    return backbone, backbone_idx
    
def order_backbone(m, backbone_idx):
    """
    Order the atoms of a molecule based on the provided backbone index and then detect the new backbone.

    :param m: The molecule object whose atoms are to be reordered.
    :type m: Chem.Mol
    :param backbone_idx: The list of indices representing the backbone atoms of the molecule.
    :type backbone_idx: list[int]
    :return: A tuple containing two elements:
        - The molecule object with its atoms reordered according to the process described.
        - The new list of indices representing the backbone atoms of the reordered molecule.
    :rtype: tuple[Chem.Mol, list[int]]

    This function first creates a copy of the provided backbone index list (`id_list`). It then iterates through all the atoms in the molecule (`m`) to find the indices of atoms that are not in the original backbone index list. These non-backbone atom indices are added to the `id_list`.

    Next, the function uses the `Chem.RenumberAtoms` method to reorder the atoms of the molecule according to the updated `id_list`. This results in a new molecule object (`m_renum`) with the atoms in a different order.

    Finally, the function calls the `detect_backbone` function (presumably defined elsewhere) on the reordered molecule (`m_renum`) to detect the new backbone. The function then returns the reordered molecule and the new list of indices representing the backbone atoms of the reordered molecule.
    """

    id_list = backbone_idx[:]
    # Renumber the atoms
    for idx in [a.GetIdx() for a in m.GetAtoms()]:
        if idx not in id_list:
            id_list.append(idx)
    m_renum = Chem.RenumberAtoms(m, newOrder=id_list)
    # Detect the backbone again
    backbone, new_backbone_idx = detect_backbone(m_renum)
    return m_renum, new_backbone_idx

def highlight_atom(m, atoms, transparency=1, isdisplay=True):
    """
    Highlight amino acid residues in cyclic peptides with different colors.

    * :py:func:`basic.detect_backbone` and :py:func:`order_backbone` function should be used before the use of cyclic peptide molecules for Identify peptide backbone and renumber atoms.

    Example::

        backbone_idx = basic.detect_backbone(m)
        m, backbone_idx = structure2sequence.order_backbone(m, backbone_idx)
        _ = structure2sequence.highlight_atom(m, [backbone_idx])

    .. image:: highlight.png
        :alt: Description of the image
        :align: center
        :scale: 30%

    """
    # Set drawing options to adjust the molecule size
    d2d = rdMolDraw2D.MolDraw2DSVG(600, 600)
    d2d.drawOptions().addAtomIndices = True
    d2d.drawOptions().addStereoAnnotation = True
    d2d.drawOptions().bondLineWidth = 2
    d2d.drawOptions().minFontSize = 15  # Adjust the font size
    d2d.drawOptions().annotationFontScale = 0.7  # Adjust the annotation font scale
    d2d.drawOptions().atomHighlightsAreCircles = True
    # Set highlight colors and radii
    highlight_colors = {}
    cmap = cm.get_cmap('tab20', len(atoms))
    for i, group in enumerate(atoms):
        for atom in group:
            if atom in highlight_colors.keys():
                highlight_colors[atom] = blend_colors(cmap(i)[:3], highlight_colors[atom], transparency)
            else:
                highlight_colors[atom] = cmap(i)[:3]  # Keep RGB, ignore alpha
    highlight_radii = {atom: 1 for group in atoms for atom in group}
    # Draw the molecule
    d2d.DrawMolecule(m, highlightAtoms=[atom for group in atoms for atom in group],
                     highlightAtomColors=highlight_colors, highlightAtomRadii=highlight_radii)
    d2d.FinishDrawing()
    # Display the image
    svg = d2d.GetDrawingText()
    # svg = re.sub(r"style='fill:([^;]+);", r"style='fill:\1;fill-opacity:"+str(transparency)+";", svg)
    if isdisplay:
        display(SVG(svg.replace('svg:', '')))
    return svg

def split_aa_unit(m):
    """
    **Identify amino acid units**

    * Use minimum amino acid unit, GLY ('NCC=O '), to recognize amino acid position on peptide;
    * Expand the branched chain to obtain complete amino acid units;
    * Amino acid units may overlap, especially at the position of ring formation;

    Example::

        aa_units = structure2sequence.split_aa_unit(m)
    
    """
    # Locate each amino acid using the minimal GLY amino acid unit
    aa_units = m.GetSubstructMatches(Chem.MolFromSmiles('NCC=O'))
    # print(aa_units)
    # Expand side chains outward from each GLY unit to obtain complete amino acid structures
    expanded_aa_set = aa_side_chain_extend(m, aa_units)
    return expanded_aa_set

# Amino acid annotation
def get_complete_aas(m, aa_units):
    """
    **Obtain complete amino acid structures**

    * ``"m"`` is generated from SMILE in mol or read in MOL file
    
    * ``"aa_units"`` is an amino acid residue separated by :py:func:`split_aa_unit` function.

    Example::

        aas, atom_mappings = struc2seq.get_complete_aas(m, aa_units)
        aa3 = struc2seq.highlight_atom(aas[3], [[]])
        IO_set.save_as_Image(aa3, 'aa3.png', 300, 300)
    
    .. image:: aa3.png
        :alt: Description of the image
        :align: center
        :scale: 30%

    """
    aas = []
    atom_mappings = []
    for atom_idxs in aa_units:
        # Create a new writable molecule object
        new_mol = RWMol()
        # Add atoms and build the index mapping
        atom_mapping = {}
        for idx in atom_idxs:
            atom = m.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        atom_mappings.append(atom_mapping)
        # Add bonds and handle external connections
        for bond in m.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_idxs and end_idx in atom_idxs:
                # Both atoms are in the list: add the bond
                new_mol.AddBond(atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType())
            elif begin_idx in atom_idxs or end_idx in atom_idxs:  # One atom is in the list, the other is not
                idx = begin_idx if begin_idx in atom_idxs else end_idx
                atom = m.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'C':
                    # Check whether the carbon connects a double-bonded oxygen and a single-bonded nitrogen (peptide bond)
                    has_double_bonded_oxygen = False
                    has_single_bonded_nitrogen = False
                    for neighbor in atom.GetNeighbors():
                        # Check whether the neighbor is an oxygen connected by a double bond
                        if neighbor.GetSymbol() == 'O' and m.GetBondBetweenAtoms(idx,
                                                                                 neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                            has_double_bonded_oxygen = True
                        elif neighbor.GetSymbol() == 'N' and m.GetBondBetweenAtoms(idx,
                                                                                   neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                            has_single_bonded_nitrogen = True
                    if has_double_bonded_oxygen and has_single_bonded_nitrogen:
                        # Add an -OH group to the carbon atom
                        new_oxygen_idx = new_mol.AddAtom(Chem.Atom('O'))
                        new_mol.AddBond(atom_mapping[idx], new_oxygen_idx, Chem.BondType.SINGLE)
                        # new_hydrogen_idx = new_mol.AddAtom(Chem.Atom('H'))
                        # new_mol.AddBond(new_oxygen_idx, new_hydrogen_idx, Chem.BondType.SINGLE)
                else:
                    pass
                    # If one atom is in the list and the other is not, add a hydrogen atom
                    # inner_idx = begin_idx if begin_idx in atom_idxs else end_idx
                    # new_idx = new_mol.AddAtom(Chem.Atom(1))  # Add a hydrogen atom as a placeholder
                    # new_mol.AddBond(atom_mapping[inner_idx], new_idx, bond.GetBondType())
        new_mol = new_mol.GetMol()
        Chem.SanitizeMol(new_mol)
        aas.append(new_mol)
    return aas, atom_mappings


def reference_aa_monomer(monomers_path):
    """
    Identify amino acids based on the monomer reference library

    * Identify the type of each amino acid unit;
    * If it can be directly recognized by 20 essential amino acids, then the amino acid can be directly determined;
    * If it cannot be directly recognized by essential amino acids, use library to find the maximum matching reference unit;
    * For cross amino acid units, select the two reference units with the most matching atoms and no cross;
    * If there are multiple best matched amino acid units, it may lead to different peptide chains;
    * Please note that some amino acids cannot be fully matched, users can view detailed information from the results;

    :param monomers_path: The "monomers_path" to the monomer library file.
    :type monomers_path: **tsv**
    :return: essential = {'code': [aa, symbol, num_atoms]}, others = {'code': [aa, symbol, num_atoms]}, 
    :rtype: dict[str, list[str]], dict[str, list[str]]

    """
    if monomers_path:
        monomers = pd.read_csv(monomers_path, sep='\t', index_col=0)
    else:
        monomers = pd.read_csv(monomer_path, sep='\t', index_col=0)

    essentials = {}
    others = {}
    temp = monomers.loc[monomers['Essential amino acids'] == 1, :]
    for i in temp.index:
        smile = temp.loc[i, 'Smiles']
        code = temp.loc[i, 'Code']
        weight = float(temp.loc[i, 'Weight'])
        symbol = temp.loc[i, 'Symbol']
        symbol = code if symbol.strip() == '' or str(symbol).upper().strip() == 'NAN' else symbol
        aa = Chem.MolFromSmiles(smile)
        num_atoms = aa.GetNumAtoms()
        essentials[code] = [aa, symbol, num_atoms]
    temp = monomers.loc[(monomers['Essential amino acids'] != 1) & monomers['Error'] != 1, :]
    for i in temp.index:
        smile = temp.loc[i, 'Smiles']
        code = temp.loc[i, 'Code']
        weight = float(temp.loc[i, 'Weight'])
        symbol = temp.loc[i, 'Symbol']
        symbol = code if str(symbol).strip() == '' or str(symbol).upper().strip() == 'NAN' else symbol
        aa = Chem.MolFromSmiles(smile)
        num_atoms = aa.GetNumAtoms()
        others[code] = [aa, symbol, num_atoms]
    return essentials, others


def get_connected_pairs(m, atom_mappings):
    """
    Verify that the main chain is connected end to end

    :param m: target molecule MOL file.
    :param atom_mappings: Cyclic peptide atom mapping.
    :type atom_mappings: obtained from :py:func:`get_complete_aas`
    :return: Amino acid-linked pair
    :rtype: List[Tuple[int, int, str]]

    """
    #print('atom', atom_mappings)
    aa_idxs = [i.keys() for i in atom_mappings]
    overlap_aa_pairs = []
    # for i in range(len(aa_idxs)-1):
    #    for j in range(i+1, len(aa_idxs)):
    #        if len(set(aa_idxs[i])&set(aa_idxs[j])):
    #            overlap_aa_pairs.append((i, j, 'side chain'))
    # Check whether the main chain is connected end to end
    for aa_i in range(len(aa_idxs) - 1):
        for aa_j in range(aa_i + 1, len(aa_idxs)):
            first_aa_idx = aa_idxs[aa_i]
            last_aa_idx = aa_idxs[aa_j]
            is_link = False
            is_peptide_bond = False
            link_count = 0
            for bond in m.GetBonds():
                begin_idx, end_idx = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
                has_peptide_bone = False
                is_side_chain = False
                if (begin_idx in first_aa_idx and end_idx in last_aa_idx) or (
                        end_idx in first_aa_idx and begin_idx in last_aa_idx):
                    is_link = True
                    link_count += 1
                    if m.GetAtomWithIdx(begin_idx).GetSymbol() == 'N' and m.GetAtomWithIdx(end_idx).GetSymbol() == 'C':
                        has_peptide_bone = False
                        is_side_chain = False
                        for neighbor in m.GetAtomWithIdx(end_idx).GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and m.GetBondBetweenAtoms(end_idx,
                                                                                     neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                has_peptide_bone = True
                            elif neighbor.GetIdx() != begin_idx and neighbor.GetIdx() in first_aa_idx and neighbor.GetIdx() in last_aa_idx:
                                is_side_chain = True
                        for neighbor in m.GetAtomWithIdx(begin_idx).GetNeighbors():
                            if neighbor.GetIdx() != end_idx and neighbor.GetIdx() in first_aa_idx and neighbor.GetIdx() in last_aa_idx:
                                is_side_chain = True

                    elif m.GetAtomWithIdx(end_idx).GetSymbol() == 'N' and m.GetAtomWithIdx(
                            begin_idx).GetSymbol() == 'C':
                        has_peptide_bone = False
                        is_side_chain = False
                        for neighbor in m.GetAtomWithIdx(begin_idx).GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and m.GetBondBetweenAtoms(begin_idx,
                                                                                     neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                has_peptide_bone = True
                            elif neighbor.GetIdx() != end_idx and neighbor.GetIdx() in first_aa_idx and neighbor.GetIdx() in last_aa_idx:
                                is_side_chain = True
                        for neighbor in m.GetAtomWithIdx(begin_idx).GetNeighbors():
                            if neighbor.GetIdx() != begin_idx and neighbor.GetIdx() in first_aa_idx and neighbor.GetIdx() in last_aa_idx:
                                is_side_chain = True
                if has_peptide_bone == True and is_side_chain == False:
                    is_peptide_bond = True
            '''if aa_i<10 and aa_j<10 and is_link:
                print('bone', aa_i, aa_j,is_peptide_bond,is_side_chain, link_count)'''
            if is_link and is_peptide_bond: 
                overlap_aa_pairs.append((aa_i, aa_j, 'peptide bond'))
                if link_count >1:
                    overlap_aa_pairs.append((aa_i, aa_j, 'side chain'))
            elif is_link: 
                overlap_aa_pairs.append((aa_i, aa_j, 'side chain'))
    #print(overlap_aa_pairs)
    return overlap_aa_pairs


def aa_matching(query_aas, atom_mappings, connected_pairs, ref_essentials, ref_others):
    """

    :param query_aas: The list of amino acids contained in the target cyclic peptide..
    :type query_aas: **List**, obtained from :py:func:`get_complete_aas`
    :param atom_mappings: Cyclic peptide atom mapping.
    :type atom_mappings: obtained from :py:func:`get_complete_aas`
    :param connected_pairs: Amino acids in cyclic peptides are linked in pairs.
    :type connected_pairs: **List**, obtained from :py:func:`get_connected_pairs`
    :param ref_essentials: The essential amino acid contained in target cyclic peptide.
    :type ref_essentials: **List**, obtained from :py:func:`reference_aa_monomer`
    :param ref_others: The nonessential amino acid contained in the target cyclic peptide.
    :type ref_others: **List**, obtained from :py:func:`reference_aa_monomer`
    :return: A treated cyclic peptide chain

    """
    result = []
    # Match against the amino acid reference library
    for aa_i in range(len(query_aas)):
        aa = query_aas[aa_i]
        aa_num_atoms = aa.GetNumAtoms()
        for ref_code, [ref, ref_symbol, ref_num_atoms] in ref_essentials.items():
            if aa_num_atoms != ref_num_atoms:
                continue
            match_idx = aa.GetSubstructMatches(ref)
            if match_idx and len(match_idx[0]) == aa_num_atoms:
                # Result: [code, query AA atom count, is essential, matched atom count, matched atom idxs, matched atom idxs in the original molecule, exact match flag]
                result.append([[ref_code, aa_num_atoms, True, len(match_idx[0]), match_idx[0], [], True]])
                break
        else:
            refs = {**ref_essentials, **ref_others}
            refs = {key: value for key, value in refs.items() if value[2] <= aa_num_atoms}
            match_codes = []
            for ref_code, [ref, ref_symbol, ref_num_atoms] in refs.items():
                match_idx = aa.GetSubstructMatches(ref)
                if match_idx:
                    match_idx = sorted(match_idx, key=lambda x: len(x), reverse=True)[0]
                    match_codes.append([ref_code, aa_num_atoms, False, len(match_idx), match_idx, [],
                                        True if aa_num_atoms == len(match_idx) else False, False])
            result.append(match_codes)
    # Remove stereoisomer duplicates, e.g., when both Cys and D-Cys match, drop D-Cys and keep Cys; also determine the atom mappings
    for res_i in range(len(result)):
        atom_mapping = atom_mappings[res_i]
        atom_mapping = {value: key for key, value in atom_mapping.items()}
        res = result[res_i]
        new_res = []
        for match_i in res:
            for match_j in res:
                if match_i[0] != match_j[0] and match_j[0] in match_i[0] and match_j[3] == match_i[3]:
                    break
            else:
                match_idx = tuple([atom_mapping[k] for k in match_i[4] if k in atom_mapping.keys()])
                match_i[5] = match_idx
                new_res.append(match_i)
        result[res_i] = new_res
    # Conflict resolution: for overlapping amino acids, find the longest non-overlapping combination
    # 1. Identify overlapping amino acids
    overlap_aa_idx = set()
    for i, j, t in connected_pairs:
        if t == 'side chain':
            overlap_aa_idx.update([i, j])
    # 2. For non-overlapping amino acids pick the largest match; for overlapping ones require non-overlapping matches covering the most atoms
    chain_aas = {}
    for aa_i in range(len(result)):
        if aa_i not in overlap_aa_idx:
            max_match_atom = 0
            max_matchs = []
            for match in result[aa_i]:
                if match[3] > max_match_atom:
                    max_match_atom = match[3]
                    max_matchs = [match]
                elif match[3] == max_match_atom:
                    max_matchs.append(match)
            chain_aas[(aa_i,)] = max_matchs
    for aa_i, aa_j, t in connected_pairs:
        if t == 'peptide bond':
            continue
        pairs = []
        for match_i in result[aa_i]:
            for match_j in result[aa_j]:
                if len(set(match_i[5]) & set(match_j[5])) == 0:
                    pairs.append([match_i, match_j, len(list(match_i[4]) + list(match_j[4]))])
        if pairs:
            max_atoms = max([i[2] for i in pairs])
        else:
            max_atoms = 0
        pairs = [i for i in pairs if i[2] >= max_atoms]
        chain_aas[(aa_i, aa_j)] = pairs
    # 3. Combine matches into chains
    chain_aas = list(chain_aas.items())
    combination_idx = [list(range(len(i[1]))) for i in chain_aas]
    combinations = list(product(*combination_idx))
    chains = []
    for combination in combinations:
        chain = [[] for i in range(len(query_aas))]
        for i in range(len(chain_aas)):
            idx = chain_aas[i][0]
            values = chain_aas[i][1]
            if len(idx) == 1:
                chain[idx[0]] = values[combination[i]]
            else:
                chain[idx[0]] = values[combination[i]][0]
                chain[idx[1]] = values[combination[i]][1]
        chains.append(chain)
    return chains


def sequence_str(chain, connected_pairs):
    """
    Generate a sequence string based on a given chain and connected pairs.

    :param chain: A list of tuples where each tuple likely represents an element in a sequence, and the first element of each tuple is used to build the sequence string.
    :type chain: list[tuple]
    :param connected_pairs: A list of tuples where each tuple contains three elements: two indices (i and j) and a connection type (t). The connection type can be 'side chain' or other types, and the indices are used to determine if a special formatting is needed for the corresponding elements in the sequence.
    :type connected_pairs: list[tuple[int, int, str]]
    :return: A string representing the sequence where elements from the chain are joined with '--' and certain elements are formatted with parentheses and a link index based on the connected pairs.
    :rtype: str

    This function first extracts the first elements of each tuple in the `chain` list to create an initial sequence list (`seq`). Then, it iterates through the `connected_pairs` list. If the connection type (`t`) is 'side chain' or the absolute difference between the two indices (`i` and `j`) is not 1, it adds parentheses and a link index to the corresponding elements in the `seq` list. Finally, it joins all the elements in the `seq` list with '--' to form the final sequence string.
    """

    seq = [i[0] for i in chain]
    link_idx = 1
    for i, j, t in connected_pairs:
        if t == 'side chain' or abs(i - j) != 1:
            seq[i] = seq[i] + '(' + str(link_idx) + ')'
            seq[j] = seq[j] + '(' + str(link_idx) + ')'
            link_idx += 1
    return '--'.join(seq)


def sequence_map(m, aa_units, chain, connected_pairs, isdisplay=True):
    """

    :param m: target molecule MOL file.
    :param aa_units: aa_units is an amino acid residue.
    :type aa_units: separated by :py:func:`split_aa_unit` function
    :param chain: cyclic peptide chain.
    :type chain: separated by :py:func:`aa_matching` function
    :param connected_pairs: Amino acids in cyclic peptides are linked in pairs.
    :type connected_pairs: **List**, obtained from :py:func:`get_connected_pairs`
    :param isdisplay: Optional "isdisplay".
    :type isdisplay: True or False

    """

    aa_locations = []
    Chem.rdDepictor.Compute2DCoords(m)
    '''try:
        m.Compute2DCoords()  # Ensure the molecule has 2D coordinates
    except:
        m = m'''
    conformer = m.GetConformer()
    for i in aa_units:
        aa_locations.append(aa_xyz(conformer, i))
    ## plot
    plt.figure(figsize=(6, 6))
    x = [i[0] for i in aa_locations]
    y = [i[1] for i in aa_locations]
    plt.scatter(x, y, color='gray', marker='o', s=1000)
    # Add numeric labels
    if aa_locations:
        for i, (xi, yi) in enumerate(zip(x, y)):
            plt.text(xi, yi, str(i + 1), color='white', fontsize=20, ha='center', va='center')
        for i, j, t in connected_pairs:
            pairx = [aa_locations[i][0], aa_locations[j][0]]
            pairy = [aa_locations[i][1], aa_locations[j][1]]
            ls = '--' if t == 'side chain' else '-'
            plt.plot(pairx, pairy, color='gray', ls=ls)
        cmap = cm.get_cmap('tab20', len(chain))
        shift = (max(x) - min(x)) * 0.02
        for aa_i in range(len(aa_units)):
            aa_symbol = chain[aa_i][0]
            # print(chain[aa_i])
            # aa_symbol = aa_symbol + '*' if chain[aa_i][1]!=chain[aa_i][-1] else aa_symbol
            plt.annotate(aa_symbol, np.array(aa_locations[aa_i]) + shift, fontsize=18, color=cmap(aa_i), zorder=999)
        # Remove the axes
        plt.axis('off')

        # Save as SVG into an in-memory string
        f = BytesIO()
        plt.savefig(f, format="svg")
        plt.close()

        # Get the SVG image string
        f.seek(0)
        svg_data = f.read().decode('utf-8')
        # Strip the namespace prefix
        svg_data = svg_data.replace('xmlns="http://www.w3.org/2000/svg" ', '')
        if isdisplay:
            display(SVG(svg_data))
        return svg_data
    else:
        print('Image failed to draw!')
        return ''


def transform(smiles, monomers_path):
    """
    Transform and report generation through an integrated function.
    
    Example::

        from IPython.display import HTML

        html = struc2seq.transform(smiles, monomers_path='monomer.tsv')
        # show html report
        HTML(html)

    :param smiles: target molecule SMILE.
    :param monomers_path: The monomers_path" to the monomer library file.
    :type monomers_path: **tsv**
    :return: transform report
    :rtype: html

    """
    report = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
        <meta charset="UTF-8">
        <title>Structure 2 Sequence Report</title>
        <style>
          .aacontainer {
            display: flex;
            flex-wrap: wrap;
            justify-content: space-around;
          }

          .aabox {
            width: 20%; /* Each box takes up about 20% of the width */
            margin: 10px;
            box-shadow: 0px 0px 10px 0px rgba(0,0,0,0.2); /* Shadow effect */
            text-align: center; /* Center inner elements */
            flex-direction: column; /* Stack children vertically */
            justify-content: center; /* Center children along the main axis */
            align-items: center; /* Center children along the cross axis */
          }

          .aabox p {
            margin: 10px 0;
            text-align: center; /* Center text */
          }

          .aabox svg {
            width: 100%;
            height: auto;
          }

          .resbox {
            width: 45%;
            margin: 10px;
            box-shadow: 0px 0px 10px 0px rgba(0,0,0,0.2); /* Shadow effect */
            text-align: center; /* Center inner elements */
            flex-direction: column; /* Stack children vertically */
            justify-content: center; /* Center children along the main axis */
            align-items: center; /* Center children along the cross axis */
          }
          .resbox p {
            margin: 10px 0;
            text-align: center; /* Center text */
          }

          .resbox svg {
            width: 100%;
            height: auto;
          }
        </style>
        </head>
        <body>
        <h1>Structure 2 Sequence Report</h1>
        <p>Structure-to-Sequence (s2s) is a computing process based on <a href='http://www.rdkit.org/'>RDkit</a> and the characteristics of cyclic peptide sequences, which can convert cyclic peptide SMILES into sequence information. This process mainly relies on the completeness of the <a href='https://www.biosino.org/iMAC/cyclicpepedia/stru2seq'>monomer reference library</a>. You can access our default monomer reference library through <a href='https://www.biosino.org/iMAC/cyclicpepedia/download'>download link</a>. The details of s2s are available on <a href='https://github.com/dfwlab/cyclicpepedia'>dfwlab/cyclicpepedia</a> on Github. And you can use this tool online on the <a href='https://www.biosino.org/iMAC/cyclicpepedia/stru2seq'>cyclicpepedia</a>.</p>
        <br/>
        <p><b>Version</b> : 1.0.1 (2023-12-26)</p>
        <hr/>
        """
    report_footer = '''
        </body>
        </html>
        '''
    # Load the structure
    sub_report = '''<h3>Load SMILES : </h3><p>SMILES : {smiles}</p><p>{loadstate}</p><hr/>'''
    try:
        m = Chem.MolFromSmiles(smiles)
        sub_report = sub_report.format(smiles=smiles, loadstate='SMILES is corrected!')
        report += sub_report
    except:
        sub_report = sub_report.format(smiles=smiles, loadstate='Load SMILES error! Check you input!')
        return report + sub_report + report_footer

    # Detect the backbone and renumber the backbone atoms
    sub_report = '''<h3>Identify peptide skeleton and renumber atoms</h3><div>{backbone}</div><hr/>'''
    is_backbone = False
    try:
        backbone, backbone_idx = detect_backbone(m)
        m, backbone_idx = order_backbone(m, backbone_idx)
        svg = highlight_atom(m, [backbone_idx], isdisplay=False)
        sub_report = sub_report.format(backbone=svg)
        report += sub_report
        is_backbone = True
    except:
        sub_report = sub_report.format(
            backbone='No single main skeleton found in the peptide. Custom amino acid sorting strategy will be used!')
        report += sub_report
        is_backbone = False
        # return report + sub_report + report_footer

    # Identify amino acid units
    sub_report = '''<h3>Identify amino acid units</h3><div>{units}</div><hr/>'''
    try:
        aa_units = split_aa_unit(m)
        svg = highlight_atom(m, aa_units, transparency=0.5, isdisplay=False)  # Overlapping parts are drawn with colors blended at fixed transparency
        sub_report = sub_report.format(units=svg)
        report += sub_report
    except:
        sub_report = sub_report.format(units='Can not find amino acid unit in peptide!')
        return report + sub_report + report_footer

    # Obtain complete amino acid structures (splitting at peptide bonds)
    aas_sub_report = '''<h3>Obtain complete amino acid structures</h3><div>{aas}</div><hr/>'''
    try:
        aas, atom_mappings = get_complete_aas(m, aa_units)
        ################ Skip outputting aas results for now ################
        # report += sub_report
    except:
        aas_sub_report = aas_sub_report.format(aas='Error!')
        return report + aas_sub_report + report_footer

    # Load the monomer reference library, identify amino acid types, and assemble amino acid chains from the molecular structure
    sub_report = '''<h3>Identify amino acids based on the monomer reference library</h3>'''
    # try:
    if monomers_path:
        ref_essentials, ref_others = reference_aa_monomer(monomers_path)
    else:
        ref_essentials, ref_others = reference_aa_monomer(monomer_path)
    connected_pairs = get_connected_pairs(m, atom_mappings)
    chains = aa_matching(aas, atom_mappings, connected_pairs, ref_essentials, ref_others)
    sub_report += '''<p>Number of chain(s) identified from the structure: <b>{nchains}</b></p>'''.format(
        nchains=len(chains))
    if not is_backbone:
        order_chain = get_connected_chain(connected_pairs, len(aas))
        aas, aa_units, chains, connected_pairs = reorder_result(order_chain, aas, aa_units, chains, connected_pairs)

    ######## Regenerate aas results (with the new numbering) ########
    temp = '<div class="aacontainer">'
    i = 1
    for aa in aas:
        svg = highlight_atom(aa, [[]], isdisplay=False)
        temp += '<div class="aabox"><p>Amino acid ' + str(i) + '</p>' + svg + '</div>'
        i += 1
    temp += '</div>'
    aas_sub_report = aas_sub_report.format(aas=temp)
    report += aas_sub_report
    ##############################################
    ci = 1
    for chain in chains:
        sub_report += '''<h4> > Chain {ci} :</h4>'''.format(ci=ci)
        sub_report += '''<p><b>Amino acid sequence :</b> {seq}</p>'''.format(seq=sequence_str(chain, connected_pairs))
        sub_report += '<div class="aacontainer">'
        sub_report += '''<div class="resbox"><p><b>Amino acid mapping</b></p>{svg}</div>'''.format(
            svg=highlight_atom(m, [i[5] for i in chain], isdisplay=False))
        sub_report += '''<div class="resbox"><p><b>Amino acid location</b></p>{svg}</div>'''.format(
            svg=sequence_map(m, aa_units, chain, connected_pairs, isdisplay=False))
        sub_report += '</div>'

        ######## Generate mapping results ########
        aas_map_sub_report = '''<p><b>Matched amino acid from monomer reference library</b></p><div>{aas}</div><hr/>'''
        refs = {**ref_essentials, **ref_others}
        temp = '<div class="aacontainer">'
        i = 1
        for aa in chain:
            ref_aa = refs[aa[0]]
            # print(aa, ref_aa)
            svg = highlight_atom(ref_aa[0], [[]], isdisplay=False)
            temp += '<div class="aabox"><p>Amino acid ' + str(i) + ': ' + aa[0] + '</p>' + svg + '</div>'
            i += 1
        temp += '</div>'
        aas_map_sub_report = aas_map_sub_report.format(aas=temp)
        sub_report += aas_map_sub_report
        ##############################################
        ci += 1
    report += sub_report
    # except:
    #    sub_report += '''<p>Error!</p>'''
    #    return report + sub_report + report_footer
    # Write the HTML content to a file
    output_file = "output.html"
    with open(output_file, "w") as f:
        f.write(report + report_footer)
    return report + report_footer


def side_chain_neighbor(m, atoms, origin_aa_idx):
    """
    Find the neighboring atoms of side chains within a molecule that are not part of the original amino acid.

    :param m: The molecule object in which to search for side chain neighbors.
    :type m: Chem.Mol
    :param atoms: A list of atom indices representing the atoms of the side chains.
    :type atoms: list[int]
    :param origin_aa_idx: A list of atom indices representing the atoms of the original amino acid.
    :type origin_aa_idx: list[int]
    :return: A list of atom indices representing the neighboring atoms of the side chains that are not part of the original amino acid.
    :rtype: list[int]

    This function first initializes an empty set (`neighbor_atoms`) to store the neighboring atom indices. It then iterates through each atom index in the `atoms` list, retrieves the corresponding atom object from the molecule (`m`), and further iterates through its neighboring atoms. For each neighboring atom, if its index is not in the `origin_aa_idx` list (meaning it's not part of the original amino acid), its index is added to the `neighbor_atoms` set. Finally, the function returns a list of atom indices from the `neighbor_atoms` set that are also not in the `atoms` list, effectively filtering out any self-neighbors and returning only the relevant neighboring atoms of the side chains that are outside the original amino acid.
    """

    neighbor_atoms = set()
    for atom_idx in atoms:
        atom = m.GetAtomWithIdx(atom_idx)
        # Iterate over this atom's neighbors
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in origin_aa_idx:
                neighbor_atoms.add(neighbor_idx)
    return [i for i in neighbor_atoms if i not in atoms]


def aa_side_chain_extend(m, aa_set):
    """
    Expand the side chains of amino acids within a molecule to include neighboring atoms.

    :param m: The molecule object containing the amino acids whose side chains are to be extended.
    :type m: Chem.Mol
    :param aa_set: A list of sets, where each set likely represents the atom indices of an amino acid.
    :type aa_set: list[set[int]]
    :return: A list of lists, where each list contains the expanded set of atom indices for an amino acid after including its side chain neighbors.
    :rtype: list[list[int]]

    This function aims to extend the side chains of amino acids in the given molecule (`m`). 

    """

    # Overlaps can occur, mainly where side chains close rings: both amino acids extend over the ring-forming side chain
    origin_aa_idx = []
    for i in aa_set:
        origin_aa_idx.extend(i[:])
    expanded_aa_set = []
    for aa in [set(i[:]) for i in aa_set]:
        neighbor_atoms = side_chain_neighbor(m, aa, origin_aa_idx)
        aa.update(neighbor_atoms)
        while (len(neighbor_atoms)):
            neighbor_atoms = side_chain_neighbor(m, aa, origin_aa_idx)
            aa.update(neighbor_atoms)
        expanded_aa_set.append(list(aa))
    return expanded_aa_set


def search_one_chain(query_idx, search_idxs, N, connected_pairs):
    """
    Search for a single chain within a set of connected pairs starting from a given query index.

    :param query_idx: The starting index for the search.
    :type query_idx: int
    :param search_idxs: A list of indices within which the search should be conducted.
    :type search_idxs: list[int]
    :param N: The maximum number of iterations or steps for the search.
    :type N: int
    :param connected_pairs: A list of tuples where each tuple contains three elements: two indices (i and j) and a connection type (t). The connection type can be used to filter out certain connections during the search.
    :type connected_pairs: list[tuple[int, int, str]]
    :return: A list representing the chain of indices found during the search, starting from the query index and following the connections specified in the connected pairs, excluding connections of type 'side chain'.
    :rtype: list[int]

    This function conducts a search to find a chain of indices within the given `search_idxs` starting from the `query_idx`. 
    """

    r = 0
    chain = [query_idx]
    while r < N:
        for i, j, t in connected_pairs:
            if i == query_idx and j in search_idxs and t != 'side chain':
                query_idx = j
                chain.append(query_idx)
                search_idxs = [k for k in search_idxs if k not in chain]
                break
            elif j == query_idx and i in search_idxs and t != 'side chain':
                query_idx = i
                chain.append(query_idx)
                search_idxs = [k for k in search_idxs if k not in chain]
                break
        r += 1
    return chain


def get_connected_chain(connected_pairs, N_aa):
    """
    Retrieve the connected chain from a set of connected pairs within a given number of amino acids.

    :param connected_pairs: A list of tuples where each tuple contains three elements: two indices (i and j) and a connection type (t). These pairs represent the connections between different elements (presumably amino acids).
    :type connected_pairs: list[tuple[int, int, str]]
    :param N_aa: The total number of amino acids in the context.
    :type N_aa: int
    :return: A list representing the connected chain of indices that covers as many of the amino acids as possible, based on the given connected pairs.
    :rtype: list[int]

    This function aims to find the connected chain within the given `connected_pairs` for a specific number of amino acids (`N_aa`). 

    """

    search_res = []
    for query_idx in range(N_aa):
        search_idxs = [i for i in range(N_aa) if i != query_idx]
        chain = search_one_chain(query_idx, search_idxs, len(search_idxs), connected_pairs)
        search_res.append(chain)
        chain = sorted(search_res, key=lambda x: len(x), reverse=True)[0]
    if len(chain) == N_aa:
        return chain
    cn = 0
    while cn < 5:  # At most 5 separate chains
        if len(chain) > N_aa - 2:
            break
        search_res = []
        for query_idx in [k for k in range(N_aa) if k not in chain]:
            search_idxs = [i for i in range(N_aa) if i != query_idx and i not in chain]
            temp = search_one_chain(query_idx, search_idxs, len(search_idxs), connected_pairs)
            search_res.append(temp)
        chain += sorted(search_res, key=lambda x: len(x), reverse=True)[0]
        cn += 1
    if len(chain) < N_aa:
        chain += [i for i in range(N_aa) if i not in chain]
    return chain


def reorder_result(order_chain, aas, aa_units, chains, connected_pairs):
    """
    Reorder various data structures related to amino acids and their connections based on a given order chain.

    :param order_chain: A list representing the new order of indices for the reordering operation.
    :type order_chain: list[int]
    :param aas: A list of amino acid related data (the specific nature of which is not clear from this function alone but is likely to be some representation of amino acids).
    :type aas: list
    :param aa_units: A list of amino acid units (again, the exact nature is not fully defined here but is related to amino acids).
    :type aa_units: list
    :param chains: A list of chains, where each chain is likely a list of indices or some other representation related to amino acids.
    :type chains: list[list]
    :param connected_pairs: A list of tuples where each tuple contains three elements: two indices (i and j) and a connection type (t), representing the connections between different elements (presumably amino acids).
    :type connected_pairs: list[tuple[int, int, str]]
    :return: A tuple containing the reordered versions of the input data structures:
        - new_aas: The reordered list of amino acid related data.
        - new_aa_units: The reordered list of amino acid units.
        - new_chains: The reordered list of chains.
        - new_connected_pairs: The reordered list of connected pairs, sorted based on the minimum of the two indices in each tuple.
    :rtype: tuple[list, list, list[list], list[tuple[int, int, str]]]

    This function takes in various data structures related to amino acids and their connections and reorders them according to the provided order chain. 

    """
    id_map = {order_chain[i]: i for i in range(len(order_chain))}
    new_aas = [None for i in range(len(aas))]
    for i in range(len(aas)):
        new_aas[id_map[i]] = aas[i]
    new_aa_units = [None for i in range(len(aa_units))]
    for i in range(len(aa_units)):
        new_aa_units[id_map[i]] = aa_units[i]
    new_chains = []
    for chain in chains:
        new_chain = [None for i in range(len(chain))]
        for j in range(len(chain)):
            new_chain[id_map[j]] = chain[j]
        new_chains.append(new_chain)
    new_connected_pairs = []
    for i, j, t in connected_pairs:
        new_i = id_map[i]
        new_j = id_map[j]
        new_connected_pairs.append((new_i, new_j, t))
    new_connected_pairs = sorted(new_connected_pairs, key=lambda x: min([x[0], x[1]]))
    return new_aas, new_aa_units, new_chains, new_connected_pairs


def aa_xyz(conformer, aa_idx):
    """
    Calculate the mean x and y coordinates of a set of atoms within a conformer corresponding to a specific amino acid index.

    :param conformer: The conformer object from which the atom positions will be retrieved.
    :type conformer: object (presumably a relevant conformer type in the context, e.g., related to molecular conformations)
    :param aa_idx: The index of the amino acid for which the atom coordinates are to be calculated.
    :type aa_idx: int
    :return: A tuple containing the mean x and y coordinates of the atoms within the conformer that belong to the specified amino acid.
    :rtype: tuple[float, float]

    """

    # Get the 2D coordinates of the specified atoms
    coords = [conformer.GetAtomPosition(idx) for idx in aa_idx]
    x_mean = sum(coord.x for coord in coords) / len(coords)
    y_mean = sum(coord.y for coord in coords) / len(coords)
    return x_mean, y_mean
