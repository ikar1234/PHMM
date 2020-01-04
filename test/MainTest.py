import unittest

from PHMM.src.Records import Record
from PHMM.src.main import PHMM
from PHMM.src.parse import parse


class MyTestCase(unittest.TestCase):
    path = "MA0466.1.meme"

    record = parse(path, ftype='meme')
    # with training data
    phmm = PHMM(record, alph="dna")

    # without training data
    phmm2 = PHMM()
    seqs = ["AATCAGACTAA",
            "ACGATCAGACT",
            "CCCATCAAACT",
            "ACGATCACTCT"]

    def test_something(self):
        self.assertEqual(True, True)

    def test_evaluate(self):
        lkh = MyTestCase.phmm.evaluate("AAAAAAAAAAA", log=False)
        self.assertEqual(lkh, 0)

        lkh = MyTestCase.phmm.evaluate("TATTGCACAAT", log=False)
        self.assertAlmostEqual(lkh, 0.036, places=3)

    # def test_viterbi_decoding(self):
    #     b = MyTestCase.phmm.viterbi_decoding("ATTGAAACAAT")
    #     print(b)

    # def test_baum_welch(self):
    #     msa = MyTestCase.phmm2.train(seqs=MyTestCase.seqs, method='baum_welch')
    #
    #     # the MSA should look like this
    #     true_msa = "A--ATCAGACT" \
    #                "ACGATCAGACT" \
    #                "CCCATCAAACT" \
    #                "ACGATC--ACT"
    #     self.assertEquals(msa, true_msa)

    def test_viterbi_training(self):
        ...

    def test_compare(self):
        # the two motifs are the same, but in different formats
        rec_jaspar = parse("../test/MA0466.1.jaspar",
                           ftype="jaspar")
        phmm_jaspar = PHMM(rec_jaspar)
        rec_meme = parse("../test/MA0466.1.meme",
                         ftype="meme")
        phmm_meme = PHMM(rec_meme)

        ratio = PHMM.compare(phmm_jaspar, phmm_meme, "TATTGCACAAT")

        self.assertAlmostEqual(ratio, 1, places=3)

        # the first motif is more distant to the sequence than the other

        rec_better = parse("../test/MA0466.1.jaspar", ftype="jaspar")
        phmm_better = PHMM(rec_better)
        print(phmm_better.Q,end='\n\n')
        rec_worse = parse("../test/MA0466.2.jaspar", ftype="jaspar")
        phmm_worse = PHMM(rec_worse)
        print(phmm_worse.Q)

        ratio = PHMM.compare(phmm_better, phmm_worse, "TATTGCACAAT")

        self.assertGreater(ratio, 1)


if __name__ == '__main__':
    unittest.main()
