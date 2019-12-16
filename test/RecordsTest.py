import unittest
import numpy as np

from ..src.Records import Record


class MyTestCase(unittest.TestCase):
    def test_something(self):
        rec = Record(name="HTSeq234", matrix=np.ones((3, 4)), ftype="MEME",
                     meta={'alength': 4, 'w': 11, 'nsites': 99494, 'E': 0})
        print(rec)
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
