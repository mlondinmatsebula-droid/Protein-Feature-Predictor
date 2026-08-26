from .amino_acid_analysis import count_amino_acids
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
    return "\n".join(lines)
