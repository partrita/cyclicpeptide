import unittest
from rdkit import Chem
from cyclicpeptide import Sequence2Structure as s2s


class TestSequence2Structure(unittest.TestCase):
    def test_essential_canonical(self):
        smi, mol = s2s.seq2stru("Ala-Ala-Cys-Asp", cyclic=True)
        self.assertIsNotNone(smi)
        self.assertIsNotNone(mol)
        self.assertEqual(mol.GetNumAtoms(), 24)
        # Check no explicit valence errors
        for a in mol.GetAtoms():
            self.assertLessEqual(a.GetExplicitValence(), 4)

    def test_essential_linear(self):
        smi, mol = s2s.seq2stru("APG", cyclic=False)
        self.assertIsNotNone(smi)
        self.assertIsNotNone(mol)
        self.assertIn("N", smi)

    def test_non_canonical_and_n_methylation(self):
        # NMe-Ala and other non-canonical residues
        smi, mol = s2s.seq2stru("cyclo[Ala-[NMe-Ala]-[4OH-Thr]-Pro]")
        self.assertIsNotNone(smi)
        self.assertIsNotNone(mol)
        # Ensure it has tertiary amide N from NMe-Ala
        tertiary_amide = Chem.MolFromSmarts("[CX3](=[OX1])N(C)[CX4]")
        self.assertTrue(mol.HasSubstructMatch(tertiary_amide))

    def test_legacy_no_essentialAA(self):
        smi, mol = s2s.seq2stru_no_essentialAA("Aad--NMe-Ala--4OH-Thr--3Me-Pro", cyclic=True)
        self.assertIsNotNone(smi)
        self.assertIsNotNone(mol)

    def test_capping_groups(self):
        # N-terminal acetyl [Ac] and C-terminal amidation [NH2]
        smi, mol = s2s.seq2stru("[Ac]-Ala-[NMe-Ala]-Gly-[NH2]", cyclic=False)
        self.assertIsNotNone(smi)
        self.assertIsNotNone(mol)
        # Acetyl cap check: CC(=O)N
        acetyl_pattern = Chem.MolFromSmarts("CC(=O)N[CX4]")
        self.assertTrue(mol.HasSubstructMatch(acetyl_pattern))
        # Amide cap check: C(=O)N
        amide_pattern = Chem.MolFromSmarts("NCC(=O)N")
        self.assertTrue(mol.HasSubstructMatch(amide_pattern))

    def test_lipid_tail_and_inline_smiles(self):
        # Octanoic acid lipid tail via [Oct]
        smi, mol = s2s.seq2stru("[Oct]-Ala-Gly", cyclic=False)
        self.assertIsNotNone(smi)
        self.assertIn("CCCCCCCC(=O)N", smi)

        # Inline SMILES
        smi2, mol2 = s2s.seq2stru("[smiles:CCCCCCCC(=O)O]-Ala-Gly", cyclic=False)
        self.assertIsNotNone(smi2)
        self.assertEqual(smi, smi2)


if __name__ == "__main__":
    unittest.main()
