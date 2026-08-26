import unittest
from src.motif_detector import detect_motifs

class TestMotifs(unittest.TestCase):
    def test_motif(self):
        seq = "NSS"
        motifs = detect_motifs(seq)
        self.assertIn('N-glycosylation', motifs)
        self.assertEqual(motifs['N-glycosylation'], [0])

if __name__ == '__main__':
    unittest.main()
