def score_seqs(seq1, seq2, match, mismatch, gap):
    """
    Write a function score_seqs(seq1, seq2, match, mismatch, gap) that computes the alignment score between two sequences (seq1 and seq2)
    using match, mismatch and gap penalty scores. Both sequences should have the same length; otherwise, the score should be 0.

    >>> score_seqs("THEFASTCAT", "THEFATCAT-", 1, -1, -2)
    -1
    >>> score_seqs("THEFASTCAT", "THEFATCA-T", 1, -1, -2)
    1
    >>> score_seqs("THEFA-TCAT", "THEFASTCAT", 1, -1, -2)
    7
    >>> score_seqs("THEFASTCAT", "THE", 1, -1, -2)
    0
    >>> score_seqs("THE-FASTCAT", "THE-FASTCAT", 1, -1, -2)
    10
    """
    len1, len2 = len(seq1), len(seq2)
    if len1 != len2:
        return 0
    

    score = sum(
        0 if i == "-" and j == "-"
        else gap if i == "-" or j == "-"
        else match if i == j
        else mismatch
        for i, j in zip(seq1, seq2)
    )
    return score
