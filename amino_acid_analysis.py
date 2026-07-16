def count_amino_acids(seq):
    counts = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    return counts
