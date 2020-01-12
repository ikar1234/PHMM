import unittest

from src.Records import Record


class MyTestCase(unittest.TestCase):
    def test_empty(self):
        empty_rec = Record(name="Empty record")
        self.assertEqual(empty_rec.alph, "dna")

    def test_unaligned(self):
        seqs = ["ACGTC", 'CAGAT', 'GCAT', 'ACTACATAC']
        empty_rec = Record(name="Unaligned record", seqs=seqs, aligned=False)
        self.assertEqual(empty_rec.alph, "dna")
        self.assertEqual(len(empty_rec.seqs), 4)

    def test_aligned(self):
        seqs = ["ACGTC",
                '-CGAT',
                'ACGACAT',
                'ACTACAT']
        empty_rec = Record(name="Unaligned record", seqs=seqs, aligned=True)
        self.assertEqual(empty_rec.alph, "dna")
        self.assertEqual(len(empty_rec.seqs), 4)


if __name__ == '__main__':
    unittest.main()
