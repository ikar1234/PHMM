import unittest

from PHMM.src.main import PHMM


class MyTestCase(unittest.TestCase):

    # def test_init(self):
    #     record = parse("better.txt", ftype='txt')
    #     # with training data
    #     phmm = PHMM(record, alph="dna")
    #     # without training data
    #     phmm2 = PHMM()

    # def test_viterbi(self):
    #     record = parse("worse.txt", ftype='txt')
    #     # with training data
    #     phmm = PHMM(record, alph="dna")

    # def test_compare(self):
    #     rec1 = parse('better.txt', ftype='txt')
    #     rec2 = parse('worse.txt', ftype='txt')
    #     phmm_better = PHMM(rec1, alph='dna')
    #     phmm_worse = PHMM(rec2, alph='dna')
    #
    #     # input is one of the sequences in the better alignment
    #     input = "CCGTGCCAAACGTAATAACACGGTA"
    #     ratio = PHMM.compare(phmm_better, phmm_worse, input)
    #     self.assertGreater(ratio, 1000)
    #     # compare the PHMM with itself
    #     ratio_one = PHMM.compare(phmm_better, phmm_better, input)
    #     self.assertAlmostEqual(ratio_one, 1, places=3)

    # def test_baum_welch(self):
    #     phmm2 = PHMM()
    #     seqs = ["ACGACT", "ACGTACT", "CGACTA", "CGATAC"]
    #     phmm2.baum_welch(seqs)

    def test_train(self):
        phmm2 = PHMM()
        seqs = ["ACGACT", "ACGTACT", "CGACTA", "CGATAC"]
        print(phmm2.train(seqs, method='baum-welch'))


if __name__ == '__main__':
    unittest.main()
