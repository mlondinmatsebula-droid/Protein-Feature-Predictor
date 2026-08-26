def composition(seq):
    counts = {}
    total = len(seq)
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1
    comp = {aa: (count/total)*100 for aa, count in counts.items()}
    return comp
