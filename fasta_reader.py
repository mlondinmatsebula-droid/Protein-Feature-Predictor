def load_fasta(filepath):
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
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
