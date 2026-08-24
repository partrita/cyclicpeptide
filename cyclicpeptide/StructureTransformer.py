from rdkit import Chem
from rdkit.Chem import AllChem

################## Structure Format Conversion ##################
# Automatically recognize and read SMILES, InChI, Molblock, SDFblock, PDBblock
# https://chemistry.stackexchange.com/questions/34563/pubchem-inchi-smiles-and-uniqueness
# InChI identifiers starting with InChI=1S/... are standard InChI; a standard InChI "must be identical for any arrangement of mobile hydrogens".
# Identifiers starting with InChI=1/... are non-standard InChI and include an extra layer beginning with /f (the fixed-H layer).
# A standard-InChI-derived structure differs from SMILES, while a non-standard InChI matches the SMILES structure; non-standard InChI generation code: Chem.MolToInchi(mol1, options='/FixedH')


# Generate SMILES, InChI, InChIKey, Molblock, SDFblock, PDBblock
def output_molecule(mol, pdbblock=None, conformation=None):
    """

    :param mol: mol file.
    :param pdbblock: Optional "conformation".
    :type pdbblock: pdbblock or None
    :param conformation: Optional "conformation".
    :type conformation: conformation or None
    :return: mol

    """

    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchi = Chem.MolToInchi(mol, options='/FixedH')
    inchikey = Chem.InchiToInchiKey(inchi)
    molblock = Chem.MolToMolBlock(mol, includeStereo=True)
    # A PDBblock needs spatial coordinates to preserve chirality, so converting other formats to PDB requires generating a 3D conformation first;
    if pdbblock is None and conformation is not None:
        pdbblock = Chem.MolToPDBBlock(conformation)
    return {'smiles': smiles, 'inchi': inchi, 'inchikey': inchikey, 'molblock': molblock, 'pdbblock': pdbblock}


def predict_3d_conformation(mol):
    """
    Predict the 3D structure of a molecule.

    :param mol: The input molecule for which the 3D conformation is to be predicted.
    :type mol: Chem.Mol
    :return: The molecule with predicted 3D conformation, including added hydrogens and generated 3D coordinates.
    :rtype: Chem.Mol

    This function takes an input molecule and predicts its 3D structure. First, it creates a copy of the input molecule and adds hydrogens to it using `Chem.AddHs`. Then, it generates 3D coordinates for the molecule with added hydrogens by employing `AllChem.EmbedMolecule` with the `AllChem.ETKDG()` method. Finally, the function returns the molecule with the predicted 3D conformation.

    """
    # molecule = Chem.MolFromSmiles(smiles)
    mol_3d = Chem.AddHs(mol.__copy__())  # Add hydrogens
    AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # Generate 3D coordinates
    return mol_3d


def mol_optimize(mol):
    """
    Optimizes the 3D structure of a molecule using the Universal Force Field (UFF).

    :param mol: The molecule to be optimized. This molecule should have 3D coordinates.
    :type mol: Chem.Mol
    :return: The optimized molecule with updated 3D coordinates.
    :rtype: Chem.Mol

    This function utilizes the Universal Force Field (UFF) to optimize the 3D structure of a given molecule. The input molecule must possess 3D coordinates for the optimization process to be carried out effectively. Once the optimization is complete, the function returns the molecule with its 3D coordinates updated to reflect the optimized structure.

    """

    # Energy minimization: UFF (Universal Force Field) suits small molecules; for peptide modeling and optimization, dedicated biomolecular tools (e.g., AMBER, CHARMM, or GROMACS) are more appropriate
    AllChem.UFFOptimizeMolecule(mol)  # Optimize with the UFF force field
    return mol



def mol2molblock(mol):
    """
    Convert a molecule object to a MolBlock string.

    :param mol: The molecule object to be converted.
    :type mol: Chem.Mol
    :return: The MolBlock string representation of the input molecule.
    :rtype: str

    """

    return Chem.MolToMolBlock(mol)


def mol2pdbblock(mol):
    """
    Convert a molecule object to a PDBBlock string.

    :param mol: The molecule object to be converted.
    :type mol: Chem.Mol
    :return: The PDBBlock string representation of the input molecule.
    :rtype: str

    """

    return Chem.MolToPDBBlock(mol)
