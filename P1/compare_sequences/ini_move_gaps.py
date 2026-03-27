def move_gaps(seq1, seq2):
    """
    Generate alignments by moving gaps in one sequence
    in another sequence.
    >>> move_gaps("THEFASTCAT", "THEFATCAT")
    ['THEFASTCAT', '-THEFATCAT', 'T-HEFATCAT', 'TH-EFATCAT', 'THE-FATCAT', 'THEF-ATCAT', 'THEFA-TCAT', 'THEFAT-CAT', 'THEFATC-AT', 'THEFATCA-T', 'THEFATCAT-']
    >>> move_gaps("THEFASTCAT", "AFASTCAT")
    ['THEFASTCAT', '--AFASTCAT', 'A--FASTCAT', 'AF--ASTCAT', 'AFA--STCAT', 'AFAS--TCAT', 'AFAST--CAT', 'AFASTC--AT', 'AFASTCA--T', 'AFASTCAT--']
    >>> move_gaps("THEFATCAT", "THEFASTCAT")
    ['THEFASTCAT', '-THEFATCAT', 'T-HEFATCAT', 'TH-EFATCAT', 'THE-FATCAT', 'THEF-ATCAT', 'THEFA-TCAT', 'THEFAT-CAT', 'THEFATC-AT', 'THEFATCA-T', 'THEFATCAT-']
    >>> move_gaps("AFASTCAT", "THEFASTCAT")
    ['THEFASTCAT', '--AFASTCAT', 'A--FASTCAT', 'AF--ASTCAT', 'AFA--STCAT', 'AFAS--TCAT', 'AFAST--CAT', 'AFASTC--AT', 'AFASTCA--T', 'AFASTCAT--']
    """
    assert len(seq1) != len(seq2)
    if len(seq1) > len(seq2):
        largest_alignment, moving_seq = seq1, seq2
    else:
        largest_alignment, moving_seq = seq2, seq1

    num_alignments = len(moving_seq)
    window = len(largest_alignment) - len(moving_seq)

    movements = [largest_alignment] + [moving_seq[:i] + "-"*window + moving_seq[i:] 
                                       for i in range(num_alignments + 1)]
    # range is num_alignments + 1 because we just want to fill the gaps in between with "-"
    
    return movements
