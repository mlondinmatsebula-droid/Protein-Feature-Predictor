#!/usr/bin/env python3
import os

files = {
    "README.md": """# Protein Feature Predictor

Analyzes protein sequences: composition, molecular weight, hydrophobicity, charge, isoelectric point, motifs, secondary structure.

## Usage
```bash
python main.py data/proteins.fasta
```
""",

    "LICENSE": "MIT License",

    "requirements.txt": "# No external dependencies",

    "main.py": """#!/usr/bin/env python3
import sys
import os
from src.fasta_reader import load_fasta
from src.report_generator import generate_protein_report

def main():
    if len(sys.argv) < 2:
        print("Usage: python main.py <fasta_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    sequences = load_fasta(fasta_file)
    report = generate_protein_report(sequences)
    print(report)
    os.makedirs('reports', exist_ok=True)
    with open('reports/protein_report.txt', 'w') as f:
        f.write(report)
    print("\\nReport saved to reports/protein_report.txt")

if __name__ == "__main__":
    main()
""",

    "src/__init__.py": "# src package",

    "src/fasta_reader.py": """def load_fasta(filepath):
    sequences = {}
    with open(filepath, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    return sequences

def save_fasta(sequences, filepath):
    with open(filepath, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\\n")
""",

    "src/amino_acid_analysis.py": """def count_amino_acids(seq):
    counts = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    return counts
""",

    "src/composition.py": """def composition(seq):
    counts = {}
    total = len(seq)
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    comp = {aa: (count/total)*100 for aa, count in counts.items()}
    return comp
""",

    "src/molecular_weight.py": """RESIDUE_WEIGHTS = {
    'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
    'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
    'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
    'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
}

def molecular_weight(seq):
    return sum(RESIDUE_WEIGHTS.get(aa, 0) for aa in seq)
""",

    "src/hydrophobicity.py": """HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

def hydrophobic_count(seq):
    return sum(1 for aa in seq if HYDROPHOBICITY.get(aa, 0) > 0)

def hydrophilic_count(seq):
    return sum(1 for aa in seq if HYDROPHOBICITY.get(aa, 0) < 0)

def longest_hydrophobic_region(seq, threshold=1.0):
    max_len = 0
    current_len = 0
    start = 0
    best_start = 0
    for i, aa in enumerate(seq):
        if HYDROPHOBICITY.get(aa, 0) > threshold:
            if current_len == 0:
                start = i
            current_len += 1
            if current_len > max_len:
                max_len = current_len
                best_start = start
        else:
            current_len = 0
    return best_start, max_len
""",

    "src/charge_analysis.py": """def charge_counts(seq):
    positive = 0
    negative = 0
    neutral = 0
    for aa in seq:
        if aa in 'RHK':
            positive += 1
        elif aa in 'DE':
            negative += 1
        else:
            neutral += 1
    return {'positive': positive, 'negative': negative, 'neutral': neutral}
""",

    "src/isoelectric_point.py": """PKA_VALUES = {
    'N_term': 8.6,
    'C_term': 3.6,
    'acidic': {'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.8},
    'basic': {'H': 6.5, 'K': 10.8, 'R': 12.5}
}

def calculate_charge(seq, ph):
    charge = 0.0
    
    # N-terminus (basic group)
    charge += 10**(PKA_VALUES['N_term'] - ph) / (1 + 10**(PKA_VALUES['N_term'] - ph))
    
    # C-terminus (acidic group)
    charge -= 1 / (1 + 10**(PKA_VALUES['C_term'] - ph))
    
    # Side chains
    for aa in seq:
        if aa in PKA_VALUES['acidic']:
            pka = PKA_VALUES['acidic'][aa]
            charge -= 1 / (1 + 10**(pka - ph))
        elif aa in PKA_VALUES['basic']:
            pka = PKA_VALUES['basic'][aa]
            charge += 10**(pka - ph) / (1 + 10**(pka - ph))
            
    return charge

def isoelectric_point(seq, tolerance=0.01):
    if not seq:
        return 7.0
        
    min_ph = 0.0
    max_ph = 14.0
    
    while True:
        mid_ph = (min_ph + max_ph) / 2
        charge = calculate_charge(seq, mid_ph)
        
        if abs(charge) <= tolerance or (max_ph - min_ph) < tolerance:
            return mid_ph
            
        if charge > 0:
            min_ph = mid_ph
        else:
            max_ph = mid_ph
""",

    "src/motif_detector.py": """import re

def detect_motifs(seq, motifs=None):
    if motifs is None:
        motifs = {
            'N-glycosylation': r'N[^P][ST][^P]',
            'PKA site': r'[RK][^P]{2}[ST]',
            'CK2 site': r'[ST][^P]{2}[DE]',
            'Phosphorylation': r'[ST]',
        }
    found = {}
    for name, pattern in motifs.items():
        matches = list(re.finditer(pattern, seq))
        if matches:
            found[name] = [m.start() for m in matches]
    return found
""",

    "src/secondary_structure.py": """def predict_secondary_structure(seq):
    helix = set('AELMQKRH')
    sheet = set('VIYCWFT')
    structure = []
    for aa in seq:
        if aa in helix:
            structure.append('H')
        elif aa in sheet:
            structure.append('E')
        else:
            structure.append('C')
    return ''.join(structure)
""",

    "src/statistics.py": """from .molecular_weight import molecular_weight
from .hydrophobicity import hydrophobic_count

def protein_statistics(sequences):
    if not sequences:
        return {}
    weights = [molecular_weight(seq) for seq in sequences.values()]
    lengths = [len(seq) for seq in sequences.values()]
    hydrophobics = [hydrophobic_count(seq) for seq in sequences.values()]
    return {
        'average_length': sum(lengths)/len(lengths),
        'longest_protein': max(lengths),
        'most_hydrophobic': max(hydrophobics),
        'highest_molecular_weight': max(weights),
    }
""",

    "src/report_generator.py": """from .amino_acid_analysis import count_amino_acids
from .composition import composition
from .molecular_weight import molecular_weight
from .hydrophobicity import hydrophobic_count, hydrophilic_count, longest_hydrophobic_region
from .charge_analysis import charge_counts
from .isoelectric_point import isoelectric_point
from .motif_detector import detect_motifs
from .secondary_structure import predict_secondary_structure

def generate_protein_report(sequences):
    lines = []
    lines.append("Protein Report")
    lines.append("=" * 50)
    lines.append("")
    for seq_id, seq in sequences.items():
        lines.append(f"Protein: {seq_id}")
        lines.append("-" * 40)
        lines.append(f"  Length: {len(seq)} aa")
        comp = composition(seq)
        lines.append("  Composition (%):")
        for aa, perc in sorted(comp.items(), key=lambda x: -x[1])[:10]:
            lines.append(f"    {aa}: {perc:.2f}")
        lines.append(f"  Molecular Weight: {molecular_weight(seq):.2f} Da")
        lines.append(f"  Hydrophobic residues: {hydrophobic_count(seq)}")
        lines.append(f"  Hydrophilic residues: {hydrophilic_count(seq)}")
        start, length = longest_hydrophobic_region(seq)
        lines.append(f"  Longest hydrophobic region: start={start}, length={length}")
        charges = charge_counts(seq)
        lines.append(f"  Charge: +{charges['positive']}, -{charges['negative']}, neutral={charges['neutral']}")
        lines.append(f"  Estimated pI: {isoelectric_point(seq):.2f}")
        motifs = detect_motifs(seq)
        if motifs:
            lines.append("  Detected motifs:")
            for motif, positions in motifs.items():
                lines.append(f"    {motif}: positions {positions[:5]}{'...' if len(positions)>5 else ''}")
        else:
            lines.append("  Detected motifs: None")
        ss = predict_secondary_structure(seq)
        helix = ss.count('H')
        sheet = ss.count('E')
        coil = ss.count('C')
        lines.append(f"  Secondary structure prediction (simple): H={helix}, E={sheet}, C={coil}")
        lines.append("")
    return "\\n".join(lines)
""",

    "src/utils.py": "# Utility functions (optional)\\npass\\n",

    "tests/test_weight.py": """import unittest
from src.molecular_weight import molecular_weight

class TestWeight(unittest.TestCase):
    def test_weight(self):
        seq = "AA"
        self.assertAlmostEqual(molecular_weight(seq), 142.1576, places=2)

if __name__ == '__main__':
    unittest.main()
""",

    "tests/test_motifs.py": """import unittest
from src.motif_detector import detect_motifs

class TestMotifs(unittest.TestCase):
    def test_motif(self):
        seq = "NSS"
        motifs = detect_motifs(seq)
        self.assertIn('N-glycosylation', motifs)
        self.assertEqual(motifs['N-glycosylation'], [0])

if __name__ == '__main__':
    unittest.main()
""",

    "tests/test_hydrophobicity.py": """import unittest
from src.hydrophobicity import hydrophobic_count, hydrophilic_count

class TestHydrophobicity(unittest.TestCase):
    def test_counts(self):
        seq = "ALV"
        self.assertEqual(hydrophobic_count(seq), 3)
        self.assertEqual(hydrophilic_count(seq), 0)

if __name__ == '__main__':
    unittest.main()
""",

    "data/proteins.fasta": """>ProteinA
MAIVMGRWKKLAVLIALAVVAGLSLG
>ProteinB
MKKLLLAVALLALVLSGCPTGG
>ProteinC
MFLVQKIVRAAVAVLLALVGA
""",
}

base = "protein-feature-predictor"
for path, content in files.items():
    full = os.path.join(base, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, 'w') as f:
        f.write(content)
    print(f"Created: {full}")

print("\\n✅ Project 2 generated. Run it with:\\n  cd protein-feature-predictor\\n  python main.py data/proteins.fasta")
