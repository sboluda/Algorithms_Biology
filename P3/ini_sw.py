from Bio.Align import substitution_matrices
subst_mat = substitution_matrices.load('BLOSUM62')

def sw(seq_i, seq_j, gap):
    """
    Local alignment of 2 sequences using Smith-Waterman algorithm
    with BLOSUM62 substitution matrix and a specified gap penalty.

    >>> sw('THEFASTCAT', 'THERAT', -4)
    Optimal score: 20.0
    ('THEFAS', 'THERAT')
    >>> sw('FAT', 'THEFASTCAT', -4)
    Optimal score: 11.0
    ('FAT', 'FAS')
    >>> sw('THECATISFAT', 'AFASTCAT', -4)
    Optimal score: 18.0
    ('CAT', 'CAT')
    """
    # Store the lengths of the sequences
    n_i, n_j = len(seq_i), len(seq_j)

    # Initialize scores and traceback matrices
    maxim, i_coord, j_coord = 0, 0, 0
    mat_scores = [[0 for _ in range(n_j + 1)] for _ in range(n_i + 1)] 
    mat_arrows = [["X" for _ in range(n_j + 1)] for _ in range(n_i + 1)] 
        
    # No need for edge-filling, given it has to be 0 (already established by default)
    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            diag = mat_scores[i-1][j-1] + subst_mat[seq_i[i-1], seq_j[j-1]]
            left = mat_scores[i][j-1] + gap # left (insertion)
            up = mat_scores[i-1][j] + gap  # up (deletion)
                
            if diag > left and diag > up and diag > 0:
                mat_scores[i][j] = int(diag)
                mat_arrows[i][j] = 0  # diagonal (match/mismatch)
                if diag > maxim:
                    maxim, i_coord, j_coord = diag, i, j
                        
            elif left > up and left > 0:
                mat_scores[i][j] = int(left)
                mat_arrows[i][j] = -1  # left (insertion)
                if left > maxim:
                    maxim, i_coord, j_coord = left, i, j

            elif up > 0:
                mat_scores[i][j] = int(up)
                mat_arrows[i][j] = 1  # up (deletion)
                if up > maxim:
                    maxim, i_coord, j_coord = up, i, j

            else: # RESET
                mat_scores[i][j] = 0
                mat_arrows[i][j] = "X"
    scores, arrows = mat_scores, mat_arrows

    # Print scores matrix (for debugging)
    '''
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print(header)
    for aa, row in zip(f" {seq_i}", scores):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))
    
    # Print traceback matrix (for debugging))
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print("\n" + header)
    for aa, row in zip(f" {seq_i}", arrows):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))
    '''

    # Print optimal score
    print(f"Optimal score: {maxim}")
    # Prepare for traceback
    aln_i, aln_j = [], []
    # Traceback - Starting at the found coordinates
    while scores[i_coord][j_coord] != 0 :  # fill in
        if arrows[i_coord][j_coord] == 0:     # diagonal
            i_coord -= 1
            j_coord -= 1
            aln_i.append(seq_i[i_coord])
            aln_j.append(seq_j[j_coord])
        elif arrows[i_coord][j_coord] == -1:  # left
            j_coord -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j_coord])
        else:        # arrows[i_coord][j_coord] == 1:   # up
            i_coord -= 1
            aln_i.append(seq_i[i_coord])
            aln_j.append("-")
    # Print the alignment
    return("".join(aln_i[::-1])), "".join(aln_j[::-1])