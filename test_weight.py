import unittest
from src.molecular_weight import molecular_weight

class TestWeight(unittest.TestCase):
    def test_weight(self):
        seq = "AA"
        self.assertAlmostEqual(molecular_weight(seq), 142.1576, places=2)

if __name__ == '__main__':
    unittest.main()
