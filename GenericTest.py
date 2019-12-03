import unittest
from collections import Counter


class MyTestCase(unittest.TestCase):
    def test_something(self):
        c = Counter(['a','a','b','b','b'])
        c = sorted(c.items(), key = lambda x: x[0])
        print(c)


if __name__ == '__main__':
    unittest.main()
