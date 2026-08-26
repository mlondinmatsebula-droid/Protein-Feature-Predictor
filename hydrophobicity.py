HYDROPHOBICITY = {
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
