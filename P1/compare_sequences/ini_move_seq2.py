from ini_score_seqs import score_seqs


def move_seq2(seq1, seq2):
    """
    >>> move_seq2("THEFASTCAT", "THEFATCAT")
    ['THEFASTCAT---------', 'THEFATCAT----------', '-THEFATCAT---------', '--THEFATCAT--------', '---THEFATCAT-------', '----THEFATCAT------', '-----THEFATCAT-----', '------THEFATCAT----', '-------THEFATCAT---', '--------THEFATCAT--', '---------THEFATCAT-', '----------THEFATCAT']
    >>> move_seq2("THEFASTCAT", "AFASTCAT")
    ['THEFASTCAT--------', 'AFASTCAT----------', '-AFASTCAT---------', '--AFASTCAT--------', '---AFASTCAT-------', '----AFASTCAT------', '-----AFASTCAT-----', '------AFASTCAT----', '-------AFASTCAT---', '--------AFASTCAT--', '---------AFASTCAT-', '----------AFASTCAT']
    >>> move_seq2("THEFASTCAT", "THECAT")
    ['THEFASTCAT------', 'THECAT----------', '-THECAT---------', '--THECAT--------', '---THECAT-------', '----THECAT------', '-----THECAT-----', '------THECAT----', '-------THECAT---', '--------THECAT--', '---------THECAT-', '----------THECAT']
    """
    assert len(seq1) != len(seq2)
    
    if len(seq1) > len(seq2):
        largest_alignment, moving_seq = seq1, seq2
    else:
        largest_alignment, moving_seq = seq2, seq1
    
    num_dashes = len(largest_alignment)
    full_seq = largest_alignment + "-"*len(moving_seq)
    movements = [full_seq] + ["-"*i + moving_seq + "-"*(num_dashes-i)
                              for i in range(len(full_seq) - len(moving_seq) + 1)]
    
    return movements


def print_scores(seq1, seq2, match, mismatch, gap):
    """
    Score and print the alignments of two sequences with match, mismatch and gap.
    >>> print_scores('THEFASTCAT', 'AFASTCAT', 1, -1, -2)
    THEFASTCAT--------
    AFASTCAT---------- -12
    -AFASTCAT--------- -12
    --AFASTCAT-------- 2
    ---AFASTCAT------- -15
    ----AFASTCAT------ -16
    -----AFASTCAT----- -19
    ------AFASTCAT---- -22
    -------AFASTCAT--- -27
    --------AFASTCAT-- -28
    ---------AFASTCAT- -33
    ----------AFASTCAT -36
    """
    list_combinations = move_seq2(seq1, seq2)
    main_align = list_combinations[0]
    print(main_align)

    for alignment in list_combinations[1:]:
        score = score_seqs(main_align, alignment, match, mismatch, gap)
        print(alignment + " " + str(score))
        
    return None