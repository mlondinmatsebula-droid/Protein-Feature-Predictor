import re

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
