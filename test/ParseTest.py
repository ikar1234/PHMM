import unittest

from parse import parse


# use this for Travis
# from src.parse import parse


class MyTestCase(unittest.TestCase):
    path = "MA0466.1.meme"
    ungapped = "ungapped.txt"
    plain = "plain.txt"

    def test_meta_alength(self):
        """
        Test that the alphabet length from the meta data is consistent with the dimension of the matrix.
        """
        p = MyTestCase.path
        rec = parse(p, ftype="meme")
        # length of sequence
        w = rec.meta.get("w", 0)
        l = rec.matrix.shape[1]

        self.assertEqual(w, l)

    def test_meta_w(self):
        """
        Test that the sequence length from the meta data is consistent with the dimension of the matrix
        """
        p = MyTestCase.path
        rec = parse(p, ftype="meme")
        # length of alphabet
        alength = rec.meta.get("alength", 0)
        a = rec.matrix.shape[0]

        self.assertEqual(alength, a)

    def test_plain(self):
        parse(MyTestCase.plain, ftype="txt")

    def test_ungapped(self):
        parse(MyTestCase.ungapped, ftype="txt")


if __name__ == '__main__':
    unittest.main()
