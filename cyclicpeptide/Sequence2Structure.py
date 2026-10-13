#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File Name: sequence2structure.py
Author: Dingfeng Wu
Creator: Dingfeng Wu
Date Created: 2024-01-03
Last Modified: 2023-01-05
Version: 1.0.0
License: MIT License
Description: Sequence-to-Structure (Seq2struc) is a computing process based on RDkit and the characteristics of cyclic peptide sequences, which can create cyclic peptide sequecne and convert sequence to cyclic peptide SMILES.

Copyright Information: Copyright (c) 2023 dfwlab (https://dfwlab.github.io/)

The code in this script can be used under the MIT License.
"""

import pandas as pd
from IPython.display import SVG, display
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import re
from .setting import *

def seq2stru(sequence, cyclic=True, monomers_path=None):
    """
    Unified Sequence-to-Structure function.
    Converts any cyclic or linear peptide sequence into SMILES and an RDKit Mol object.
    
    Supports:
      - Essential amino acids (1-letter e.g. 'APG' or 3-letter e.g. 'Ala-Pro-Gly')
      - Non-canonical amino acids in brackets (e.g. '[NMe-Ala]', '[4OH-Thr]', '[Aad]', '[3Me-Pro]')
      - N-methylation of peptide bonds
      - Other chemical moieties:
          - Terminal capping groups: '[Ac]' / 'Ac' (acetyl), '[NH2]' / 'NH2' (amidation), '[Oct]' (lipid tail)
          - Inline SMILES: '[smiles:CCCCCCCC(=O)O]'
      - Both cyclic (cyclic=True or 'cyclo[...]' prefix) and linear peptides.
      
    Example::
    
        smiles, mol = Sequence2Structure.seq2stru('cyclo[Ala-[NMe-Ala]-[4OH-Thr]-Pro]')
        print(smiles)
    """
    try:
        s = sequence.strip()
        is_cyclo = False
        if s.startswith("cyclo[") and s.endswith("]"):
            is_cyclo = True
            s = s[6:-1]
        elif "cyclo" in s:
            is_cyclo = True
            s = s.replace("cyclo", "").strip("[]")
        cyclic = cyclic or is_cyclo
        s = re.sub(r"\(\d+\)", "", s)
        
        tokens = []
        if "--" in s:
            tokens = [t.strip() for t in s.split("--") if t.strip()]
        elif "[" in s or "-" in s:
            s_masked = re.sub(r"\[[^\]]*\]", lambda m: m.group(0).replace("-", "~~"), s)
            for part in s_masked.split("-"):
                part = part.replace("~~", "-").strip()
                if part:
                    tokens.append(part)
        elif len(s) > 0:
            if all(c.isalpha() and c.isupper() for c in s):
                tokens = list(s)
            else:
                tokens = [s]
                
        if monomers_path:
            monomers = pd.read_csv(monomers_path, sep="\t", index_col=0)
        else:
            monomers = pd.read_csv(monomer_path, sep="\t", index_col=0)
        if "Error" in monomers.columns:
            monomers = monomers.loc[monomers["Error"] != 1, :]
            
        capping_dict = {
            "AC": "CC(=O)O",
            "NAC": "CC(=O)O",
            "ACETYL": "CC(=O)O",
            "NH2": "N",
            "AMIDE": "N",
            "OCT": "CCCCCCCC(=O)O",
            "PAM": "CCCCCCCCCCCCCCCC(=O)O",
        }
        
        one_to_three = {}
        try:
            with open(AminoAcids_path, "r", encoding="utf-8") as f:
                for line in f:
                    item = eval(line.strip())
                    one_to_three[item[0].upper()] = item[1]
        except:
            pass
            
        monomer_dict = {}
        for _, row in monomers.iterrows():
            c = str(row["Code"]).strip().upper()
            smi_str = str(row["Smiles"]).strip()
            if c and smi_str:
                monomer_dict[c] = smi_str
            n = str(row["Name"]).strip().upper()
            if n and smi_str and n not in monomer_dict:
                monomer_dict[n] = smi_str
                
        def _resolve(token):
            clean = token.strip()
            if clean.startswith("[") and clean.endswith("]"):
                clean = clean[1:-1].strip()
            if clean.lower().startswith("smiles:"):
                return Chem.MolFromSmiles(clean[7:].strip())
            up = clean.upper()
            if up in capping_dict:
                return Chem.MolFromSmiles(capping_dict[up])
            if len(clean) == 1 and up in one_to_three:
                clean = one_to_three[up]
                up = clean.upper()
            if up in monomer_dict:
                return Chem.MolFromSmiles(monomer_dict[up])
            return None
            
        mols = [_resolve(t) for t in tokens]
        if None in mols:
            return None, None
            
        p = build_peptide_from_mols(mols, cyclic=cyclic)
        if p is None:
            return None, None
        p = Chem.RemoveHs(p)
        Chem.AssignAtomChiralTagsFromStructure(p)
        smiles = Chem.MolToSmiles(p, canonical=True)
        return smiles, p
    except:
        return None, None


def seq2stru_essentialAA(sequence, cyclic=True):
    try:
        smi, mol = seq2stru(sequence, cyclic=cyclic)
        if smi is not None:
            return smi, mol
        if '-' in sequence:
            sequence = [code2symbol(code)[0] for code in sequence.split('-')]
            sequence = ''.join(sequence)
        peptide = create_peptide_of_essentialAA(sequence, cyclic=cyclic)
        peptide = Chem.RemoveHs(peptide)
        Chem.AssignAtomChiralTagsFromStructure(peptide)
        smiles = Chem.MolToSmiles(peptide, canonical=True)
        return smiles, peptide
    except:
        return None, None


def seq2stru_no_essentialAA(sequence, references=None, cyclic=True):
    """
    Example::

        sequence = 'Aad--NMe-Ala--4OH-Thr--3Me-Pro'
        smiles, peptide = sequence2structure.seq2stru_no_essentialAA(sequence, references, cyclic=True)
        print(smiles)
        IOManager.plot_smiles(smiles, w=300, h=300, isdisplay=True)

    SMILES : CC1CCN2C(=O)C(C(O)CO)NC(=O)C(C)N(C)C(=O)C(CCCC(=O)O)NC(=O)C12
    """
    try:
        smi, mol = seq2stru(sequence, cyclic=cyclic)
        if smi is not None:
            return smi, mol
        if references is None:
            references = reference_aa_monomer(None)
        seq_tokens = [i.strip() for i in re.sub(r'\(\d+\)', '', sequence).split('--')]
        peptide = create_peptide_of_no_essentialAA(seq_tokens, references, cyclic=cyclic)
        smiles = Chem.MolToSmiles(peptide, canonical=True)
        try:
            peptide = Chem.RemoveHs(peptide)
            Chem.AssignAtomChiralTagsFromStructure(peptide)
            return smiles, peptide
        except:
            return smiles, peptide
    except:
        return None, None


def code2symbol(code):
    """

    :param code: sequence.
    :type code: amino acid chain format
    :return: One letter code

    """

    # Open the text file and read its contents
    file_path = AminoAcids_path  # Replace with your file path
    AminoAcids = []
    with open(file_path, 'r', encoding='utf-8') as file:
        """Reads the contents of the file line by line and adds it to the list"""
        for line in file:
            """Strip the trailing newline and append the entry to the list"""
            AminoAcids.append(eval(line.strip()))
    c2s = {i[1]: [i[0], i[2]] for i in AminoAcids}
    if code.upper() in c2s.keys():
        return c2s[code.upper()]
    else:
        return None


def find_terminals(mol):
    """
    Find the backbone N-terminal (amine nitrogen) and C-terminal (carboxylic acid carbon).
    Returns (n_idx, c_idx).
    """
    if mol is None:
        return None, None
    c_patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    c_matches = [m[0] for m in mol.GetSubstructMatches(c_patt)]
    n_patt = Chem.MolFromSmarts("[NX3;H1,H2,H0,H3;!$(N=*);!$(N~O)]")
    n_matches = [m[0] for m in mol.GetSubstructMatches(n_patt)]
    
    if not c_matches and not n_matches:
        return None, None
    if not c_matches:
        return (n_matches[0] if n_matches else None, None)
    if not n_matches:
        return (None, c_matches[0] if c_matches else None)
        
    best_pair = None
    min_dist = 999
    for n in n_matches:
        for c in c_matches:
            path = Chem.GetShortestPath(mol, n, c)
            if path and 1 < len(path) <= 4:
                if len(path) < min_dist:
                    min_dist = len(path)
                    best_pair = (n, c)
    if best_pair:
        return best_pair
    return (n_matches[0], c_matches[0])


def build_peptide_from_mols(mol_list, cyclic=True):
    """
    Build a peptide molecule from a list of monomer Mol objects.
    Uses atom tagging to prevent valence errors and index shifting bugs.
    """
    if not mol_list:
        return None
    tagged_mols = []
    for i, m in enumerate(mol_list):
        m_copy = Chem.RWMol(m)
        n_idx, c_idx = find_terminals(m_copy)
        for a in m_copy.GetAtoms():
            a.SetIntProp("_res_id", i)
        if n_idx is not None:
            m_copy.GetAtomWithIdx(n_idx).SetIntProp("_is_n_term", 1)
        if c_idx is not None:
            m_copy.GetAtomWithIdx(c_idx).SetIntProp("_is_c_term", 1)
        tagged_mols.append(m_copy)
        
    combined = tagged_mols[0]
    for m in tagged_mols[1:]:
        combined = Chem.CombineMols(combined, m)
    combined = Chem.RWMol(combined)
    
    to_delete = set()
    bonds_to_add = []
    for i in range(len(mol_list) - 1):
        c_atom_idx = None
        n_atom_idx = None
        for a in combined.GetAtoms():
            if a.GetIntProp("_res_id") == i and a.HasProp("_is_c_term"):
                c_atom_idx = a.GetIdx()
            elif a.GetIntProp("_res_id") == i + 1 and a.HasProp("_is_n_term"):
                n_atom_idx = a.GetIdx()
                
        if c_atom_idx is not None and n_atom_idx is not None:
            c_atom = combined.GetAtomWithIdx(c_atom_idx)
            for nb in c_atom.GetNeighbors():
                if nb.GetSymbol() == "O" and combined.GetBondBetweenAtoms(c_atom_idx, nb.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                    to_delete.add(nb.GetIdx())
                    for h in nb.GetNeighbors():
                        if h.GetAtomicNum() == 1:
                            to_delete.add(h.GetIdx())
                    break
            bonds_to_add.append((i, i + 1, "_is_c_term", "_is_n_term"))
            
    if cyclic and len(mol_list) >= 2:
        last_i = len(mol_list) - 1
        first_i = 0
        c_atom_idx = None
        n_atom_idx = None
        for a in combined.GetAtoms():
            if a.GetIntProp("_res_id") == last_i and a.HasProp("_is_c_term"):
                c_atom_idx = a.GetIdx()
            elif a.GetIntProp("_res_id") == first_i and a.HasProp("_is_n_term"):
                n_atom_idx = a.GetIdx()
        if c_atom_idx is not None and n_atom_idx is not None:
            c_atom = combined.GetAtomWithIdx(c_atom_idx)
            for nb in c_atom.GetNeighbors():
                if nb.GetSymbol() == "O" and combined.GetBondBetweenAtoms(c_atom_idx, nb.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                    to_delete.add(nb.GetIdx())
                    for h in nb.GetNeighbors():
                        if h.GetAtomicNum() == 1:
                            to_delete.add(h.GetIdx())
                    break
            bonds_to_add.append((last_i, first_i, "_is_c_term", "_is_n_term"))
            
    for idx in sorted(to_delete, reverse=True):
        combined.RemoveAtom(idx)
        
    for res_from, res_to, prop_from, prop_to in bonds_to_add:
        from_idx = None
        to_idx = None
        for a in combined.GetAtoms():
            if a.GetIntProp("_res_id") == res_from and a.HasProp(prop_from):
                from_idx = a.GetIdx()
            elif a.GetIntProp("_res_id") == res_to and a.HasProp(prop_to):
                to_idx = a.GetIdx()
        if from_idx is not None and to_idx is not None:
            if not combined.GetBondBetweenAtoms(from_idx, to_idx):
                combined.AddBond(from_idx, to_idx, Chem.BondType.SINGLE)
                
    for a in combined.GetAtoms():
        for p in ["_res_id", "_is_n_term", "_is_c_term"]:
            if a.HasProp(p):
                a.ClearProp(p)
                
    res = combined.GetMol()
    try:
        Chem.SanitizeMol(res)
        Chem.AssignAtomChiralTagsFromStructure(res)
    except:
        pass
    return res


def link_aa_by_peptide_bond(mol, c_index, n_index):
    """
    Link two amino acids by a peptide bond in a given molecule.
    Uses atom tagging to safely remove leaving OH groups without index shift corruption.
    """
    m = Chem.RWMol(mol)
    m.GetAtomWithIdx(c_index).SetIntProp("_target_c", 1)
    m.GetAtomWithIdx(n_index).SetIntProp("_target_n", 1)
    
    to_delete = []
    c_atom = m.GetAtomWithIdx(c_index)
    for neighbor in c_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'O' and m.GetBondBetweenAtoms(c_index, neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
            to_delete.append(neighbor.GetIdx())
            for h_atom in neighbor.GetNeighbors():
                if h_atom.GetAtomicNum() == 1:
                    to_delete.append(h_atom.GetIdx())
            break
            
    for idx in sorted(set(to_delete), reverse=True):
        m.RemoveAtom(idx)
        
    c_final = None
    n_final = None
    for atom in m.GetAtoms():
        if atom.HasProp("_target_c"):
            c_final = atom.GetIdx()
        elif atom.HasProp("_target_n"):
            n_final = atom.GetIdx()
            
    if c_final is not None and n_final is not None:
        if not m.GetBondBetweenAtoms(c_final, n_final):
            m.AddBond(c_final, n_final, order=Chem.rdchem.BondType.SINGLE)
        m.GetAtomWithIdx(c_final).ClearProp("_target_c")
        m.GetAtomWithIdx(n_final).ClearProp("_target_n")
        
    res = m.GetMol()
    try:
        Chem.SanitizeMol(res)
    except:
        pass
    return res


def detect_backbone(mol):
    """
    Detect the backbone of a given molecule.
    """
    bbsmiles = "C(=O)CN" * len(
        mol.GetSubstructMatches(Chem.MolFromSmiles('NCC=O')))
    matches = mol.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))
    backbone = ''
    if matches:
        backbone = matches[0]
    else:
        # Fallback to peptide bond matching
        patt = Chem.MolFromSmarts("[CX3](=[OX1])-[NX3;!$(N=*)]")
        pb_matches = mol.GetSubstructMatches(patt)
        if pb_matches:
            bb_atoms = []
            for m in pb_matches:
                bb_atoms.extend([m[0], m[2]])
            backbone = tuple(dict.fromkeys(bb_atoms))
        else:
            return backbone, None
    backbone_idx = list(backbone)
    backbone_idx.reverse()
    return backbone, backbone_idx


def reference_aa_monomer(monomers_path=None): 
    """
    Load monomer library from monomer.tsv file.
    """
    if monomers_path:
        monomers = pd.read_csv(monomers_path, sep='\t', index_col=0)
    else:
        monomers = pd.read_csv(monomer_path, sep='\t', index_col=0)
    if 'Error' in monomers.columns:
        monomers = monomers.loc[monomers['Error'] != 1, :]
    references = {}
    for i in monomers.index:
        smile = monomers.loc[i, 'Smiles']
        code = str(monomers.loc[i, 'Code']).strip()
        try:
            weight = float(monomers.loc[i, 'Weight'])
        except:
            weight = 0.0
        symbol = str(monomers.loc[i, 'Symbol']).strip()
        symbol = code if symbol == '' or symbol.upper() == 'NAN' else symbol
        try:
            aa = Chem.MolFromSmiles(smile)
            if aa is None:
                continue
            n_idx, c_idx = find_terminals(aa)
            if n_idx is not None and c_idx is not None:
                num_atoms = aa.GetNumAtoms()
                references[code] = [aa, symbol, num_atoms, [n_idx, 0, c_idx, 0]]
        except:
            pass
    return references


def create_peptide_of_no_essentialAA(sequence, references=None, cyclic=True):
    """
    Create a peptide molecule from non-essential / modified amino acid sequence.
    """
    if references is None:
        references = reference_aa_monomer(None)
    mols = []
    for aa_code in sequence:
        if aa_code in references:
            mols.append(references[aa_code][0])
        else:
            m = Chem.MolFromSmiles(aa_code)
            if m is not None:
                mols.append(m)
    if mols:
        return build_peptide_from_mols(mols, cyclic=cyclic)
    return None


def connect_two_aa_with_peptide_bond(aa1, c_index1, aa2, c_index2, n_index2):
    """
    Connect two amino acids with a peptide bond using tag-based atom tracking.
    """
    m1 = Chem.RWMol(aa1)
    m2 = Chem.RWMol(aa2)
    m1.GetAtomWithIdx(c_index1).SetIntProp("_c1", 1)
    m2.GetAtomWithIdx(n_index2).SetIntProp("_n2", 1)
    m2.GetAtomWithIdx(c_index2).SetIntProp("_c2", 1)
    
    combined = Chem.RWMol(Chem.CombineMols(m1, m2))
    
    c1_idx = None
    to_delete = []
    for a in combined.GetAtoms():
        if a.HasProp("_c1"):
            c1_idx = a.GetIdx()
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == 'O' and combined.GetBondBetweenAtoms(c1_idx, nb.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                    to_delete.append(nb.GetIdx())
                    for h in nb.GetNeighbors():
                        if h.GetAtomicNum() == 1:
                            to_delete.append(h.GetIdx())
                    break
            break
            
    for idx in sorted(set(to_delete), reverse=True):
        combined.RemoveAtom(idx)
        
    c1_final = None
    n2_final = None
    c2_final = None
    for a in combined.GetAtoms():
        if a.HasProp("_c1"):
            c1_final = a.GetIdx()
        elif a.HasProp("_n2"):
            n2_final = a.GetIdx()
        if a.HasProp("_c2"):
            c2_final = a.GetIdx()
            
    if c1_final is not None and n2_final is not None:
        if not combined.GetBondBetweenAtoms(c1_final, n2_final):
            combined.AddBond(c1_final, n2_final, Chem.BondType.SINGLE)
            
    if c1_final is not None:
        combined.GetAtomWithIdx(c1_final).ClearProp("_c1")
    if n2_final is not None:
        combined.GetAtomWithIdx(n2_final).ClearProp("_n2")
    if c2_final is not None:
        combined.GetAtomWithIdx(c2_final).ClearProp("_c2")
        
    res = combined.GetMol()
    try:
        Chem.SanitizeMol(res)
    except:
        pass
    return res, c2_final


def get_molblock(mol):
    """
    Retrieve the MolBlock representation of a given molecule.

    :param mol: The molecule object for which the MolBlock is to be retrieved.
    :type mol: Chem.Mol
    :return: The MolBlock string representation of the molecule. The `forceV3000=True` parameter indicates that the MolBlock should be in the V3000 format (if available and applicable).
    :rtype: str

    """

    return Chem.MolToMolBlock(mol, forceV3000=True)


def get_molblock_from_smiles(smiles):
    """
    Generate the MolBlock representation of a molecule from its SMILES string.

    :param smiles: The SMILES string representing the molecule.
    :type smiles: str
    :return: The MolBlock string representation of the molecule. The `forceV3000=True` indicates that the MolBlock will be in the V3000 format (if available and applicable).
    :rtype: str

    This function first converts the given SMILES string into a molecule object using `Chem.MolFromSmiles`. If the conversion is successful, it then retrieves the MolBlock representation of the molecule using `Chem.MolToMolBlock` with the `forceV3000=True` option to potentially get the MolBlock in the V3000 format. If the conversion from SMILES to molecule fails, it will likely return `None` as `Chem.MolToMolBlock` will be called on a `None` object (since `mol` would be `None` if the `Chem.MolFromSmiles` conversion fails).
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToMolBlock(mol, forceV3000=True)
