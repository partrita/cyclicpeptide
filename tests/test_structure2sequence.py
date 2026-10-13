import unittest
from rdkit import Chem
from cyclicpeptide import Structure2Sequence as struc2seq
from cyclicpeptide import Sequence2Structure as s2s


class TestStructure2Sequence(unittest.TestCase):
    def test_reference_loading(self):
        ess, others = struc2seq.reference_aa_monomer(None)
        self.assertGreater(len(ess), 0)
        self.assertGreater(len(others), 0)
        # Ensure NMe-Ala and 3Me-Pro are in others
        self.assertIn("NMe-Ala", others)
        self.assertIn("3Me-Pro", others)

    def test_transform_non_canonical(self):
        smiles = "CC1CCN2C(=O)C(C(O)CO)NC(=O)C(C)N(C)C(=O)C(CCCC(=O)O)NC(=O)C12"
        html = struc2seq.transform(smiles, monomers_path=None)
        self.assertIsNotNone(html)
        self.assertIn("Structure 2 Sequence Report", html)
        # Should identify 3Me-Pro, NMe-Ala, Aad, 4OH-Thr
        self.assertIn("[3Me-Pro]", html)
        self.assertIn("[NMe-Ala]", html)
        self.assertIn("[Aad]", html)
        self.assertIn("[4OH-Thr]", html)

    def test_transform_essential(self):
        smiles, _ = s2s.seq2stru("Ala-Ala-Cys-Asp", cyclic=True)
        html = struc2seq.transform(smiles, monomers_path=None)
        self.assertIsNotNone(html)
        self.assertIn("Ala", html)
        self.assertIn("Cys", html)
        self.assertIn("Asp", html)


if __name__ == "__main__":
    unittest.main()
