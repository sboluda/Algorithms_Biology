from ini_move_gaps import move_gaps
from ini_score_seqs import score_seqs


def move_gaps_scores(seq1, seq2, match, mismatch, gap):
    """
    >>> move_gaps_scores('THEFASTCAT', 'AFASTCAT', 1, -1, -2)
    THEFASTCAT
    --AFASTCAT 2
    A--FASTCAT 2
    AF--ASTCAT 0
    AFA--STCAT -2
    AFAS--TCAT -4
    AFAST--CAT -6
    AFASTC--AT -8
    AFASTCA--T -10
    AFASTCAT-- -12
    >>> move_gaps_scores('AFASTCAT', 'THEFASTCAT', 1, -1, -2)
    THEFASTCAT
    --AFASTCAT 2
    A--FASTCAT 2
    AF--ASTCAT 0
    AFA--STCAT -2
    AFAS--TCAT -4
    AFAST--CAT -6
    AFASTC--AT -8
    AFASTCA--T -10
    AFASTCAT-- -12
    """
    list_combinations = move_gaps(seq1, seq2)
    main_align = list_combinations[0]
    print(main_align)

    for alignment in list_combinations[1:]:
        score = score_seqs(main_align, alignment, match, mismatch, gap)
        print(alignment + " " + str(score))
        
    return None
