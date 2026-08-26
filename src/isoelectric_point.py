PKA_VALUES = {
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
