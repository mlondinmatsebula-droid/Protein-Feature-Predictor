#!/usr/bin/env python3
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
    print("\nReport saved to reports/protein_report.txt")

if __name__ == "__main__":
    main()
