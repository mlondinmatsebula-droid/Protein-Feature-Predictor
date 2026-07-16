import unittest
from src.hydrophobicity import hydrophobic_count, hydrophilic_count

class TestHydrophobicity(unittest.TestCase):
    def test_counts(self):
        seq = "ALV"
        self.assertEqual(hydrophobic_count(seq), 3)
        self.assertEqual(hydrophilic_count(seq), 0)

if __name__ == '__main__':
    unittest.main()
