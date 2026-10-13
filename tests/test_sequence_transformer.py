import unittest
from cyclicpeptide import SequenceTransformer as st


class TestSequenceTransformer(unittest.TestCase):
    def test_canonical_iupac(self):
        seq = "Ala-Gly-Pro"
        fmt, nodes, edges = st.read_sequence(seq)
        self.assertEqual(fmt, "IUPAC condensed")
        self.assertEqual(nodes, ["Ala", "Gly", "Pro"])
        self.assertEqual(edges, [])

    def test_one_letter(self):
        seq = "AGP"
        fmt, nodes, edges = st.read_sequence(seq)
        self.assertEqual(fmt, "One letter peptide")
        self.assertEqual(nodes, ["ALA", "GLY", "PRO"])

    def test_amino_acid_chain(self):
        seq = "Aad--NMe-Ala--4OH-Thr--3Me-Pro"
        fmt, nodes, edges = st.read_sequence(seq)
        self.assertEqual(fmt, "Amino acid chain")
        self.assertEqual(nodes, ["Aad", "NMe-Ala", "4OH-Thr", "3Me-Pro"])

    def test_bracketed_non_canonical(self):
        seq = "Ala-[NMe-Ala]-[4OH-Thr]-Pro"
        fmt, nodes, edges = st.read_sequence(seq)
        self.assertEqual(fmt, "IUPAC condensed")
        self.assertEqual(nodes, ["Ala", "[NMe-Ala]", "[4OH-Thr]", "Pro"])

    def test_cyclo_bracketed(self):
        seq = "cyclo[Ala-[NMe-Ala]-[4OH-Thr]-Pro]"
        fmt, nodes, edges = st.read_sequence(seq)
        self.assertEqual(fmt, "IUPAC condensed")
        self.assertEqual(nodes, ["Ala", "[NMe-Ala]", "[4OH-Thr]", "Pro"])
        self.assertEqual(edges, [(0, 3)])

    def test_roundtrip_creation(self):
        nodes = ["Ala", "[NMe-Ala]", "[4OH-Thr]", "Pro"]
        edges = [(0, 3)]
        seq_dict = st.create_sequence(nodes, edges)
        self.assertEqual(seq_dict["iupac_condensed"], "cyclo[Ala-[NMe-Ala]-[4OH-Thr]-Pro]")
        self.assertIn("Ala(1)--[NMe-Ala]--[4OH-Thr]--Pro(1)", seq_dict["amino_acid_chain"])


if __name__ == "__main__":
    unittest.main()
