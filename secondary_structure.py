def predict_secondary_structure(seq):
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
