import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)

    def test_evaluate(self):
        ...

    def test_viterbi_decoding(self):
        ...

    def test_baum_welch(self):
        ...

    def test_viterbi_training(self):
        ...


if __name__ == '__main__':
    unittest.main()
