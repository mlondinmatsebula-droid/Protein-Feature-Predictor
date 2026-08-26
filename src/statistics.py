from .molecular_weight import molecular_weight
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
