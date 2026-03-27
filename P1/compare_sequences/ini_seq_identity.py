def seq_identity(seq1, seq2):
    """
    Return the percentage of sequence identity. Exclude positions with gaps on any sequence
    >>> seq_identity("FASTCAT", "FATCAT")

    >>> seq_identity("FASTCAT", "FASTCAT")
    100.0
    >>> seq_identity("FASTCAT", "FASTRAT")
    85.7
    >>> seq_identity("-FASTCAT", "-FASTRAT")
    85.7
    >>> seq_identity("FASTCAT", "FA-TCAT")
    100.0
    >>> seq_identity("FASTCAT", "FAT-CAT")
    83.3
    >>> seq_identity("AFASTCAT", "-FASTRAT")
    85.7
    >>> seq_identity("FASTCAT", "AAAAAAA")
    28.6
    >>> seq_identity("FASTCAT", "AFAAAFA")
    0.0
    """
    if len(seq1) != len(seq2): 
        return None
    pairs = list(zip(seq1, seq2))
    matches = sum(1 for x, y in pairs if "-" not in (x, y) and x == y) 
    length = sum(1 for x, y in pairs if "-" not in (x, y))
    return round((matches/length)*100, 1)