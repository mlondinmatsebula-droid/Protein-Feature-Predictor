def charge_counts(seq):
    positive = 0
    negative = 0
    neutral = 0
    for aa in seq:
        if aa in 'RHK':
            positive += 1
        elif aa in 'DE':
            negative += 1
        else:
            neutral += 1
    return {'positive': positive, 'negative': negative, 'neutral': neutral}
