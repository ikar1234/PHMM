import unittest

from PHMM.src.Records import Record
from PHMM.src.main import PHMM
from PHMM.src.parse import parse


class MyTestCase(unittest.TestCase):
    path = "MA0466.1.meme"

    record = parse(path, ftype='meme')
    # with training data
    # phmm = PHMM(record, alph="dna")

    # without training data
    phmm2 = PHMM()
    seqs = ["AATAGACTAA",
            "ACGATCAGACT",
            "CCCATCAAACT",
            "ACGATCACTCT"]

    # def test_something(self):
    #     rec = parse('plain.txt',ftype='txt')
    #     pl = PHMM(rec)
    #     print(pl.P)
    #     print(pl.Q)

    def test_make_paths(self):
        # Examlpe is taken from Durbin
        states = [0, 0, 1, 1, 0]
        seqs = ["VG--H", "V---N", "VE--D", "IAADN"]

        p = PHMM.make_paths(states, seqs)
        print(list(p))

    # def test_evaluate(self):
    #     lkh = MyTestCase.phmm.evaluate("AAAAAAAAAAA", log=False)
    #     self.assertEqual(lkh, 0)
    #
    #     lkh = MyTestCase.phmm.evaluate("TATTGCACAAT", log=False)
    #     self.assertAlmostEqual(lkh, 0.036, places=3)
    #
    # def test_viterbi_decoding(self):
    #     gapped_phmm = PHMM(parse(filepath='../test/plain.txt', ftype='txt'))
    #     b = list(gapped_phmm.viterbi_decoding("ATTGACGTACGTAAT"))
    #     print(b)
    #
    # def test_baum_welch(self):
    #     gapped_phmm = PHMM(parse(filepath='../test/plain.txt', ftype='txt'))
    #
    #     msa = gapped_phmm.train(seqs=MyTestCase.seqs, method='viterbi')
    #     print(msa)

    # def test_compare(self):
    #     seq = "TATTGCACAAT"
    #     # the two motifs are the same, but in different formats
    #     rec_jaspar = parse("../test/MA0466.1.jaspar",
    #                        ftype="jaspar")
    #     phmm_jaspar = PHMM(rec_jaspar)
    #     rec_meme = parse("../test/MA0466.1.meme",
    #                      ftype="meme")
    #     phmm_meme = PHMM(rec_meme)
    #
    #     ratio = PHMM.compare(phmm_jaspar, phmm_meme, seq)
    #
    #     self.assertAlmostEqual(ratio, 1, places=3)
    #
    #     # the first motif is more distant to the sequence than the other
    #
    #     rec_better = parse("../test/MA0466.1.jaspar", ftype="jaspar")
    #     phmm_better = PHMM(rec_better)
    #     print(phmm_better._forward(seq), end='\n\n')
    #     rec_worse = parse("../test/MA0466.2.jaspar", ftype="jaspar")
    #     phmm_worse = PHMM(rec_worse)
    #     print(phmm_worse._forward(seq))
    #
    #     ratio = PHMM.compare(phmm_better, phmm_worse, seq)
    #
    #     self.assertGreater(ratio, 1)


if __name__ == '__main__':
    unittest.main()
