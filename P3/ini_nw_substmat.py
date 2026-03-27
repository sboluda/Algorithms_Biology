from Bio.Align import substitution_matrices
subst_mat = substitution_matrices.load('BLOSUM62')

def nw(seq_i, seq_j, gap):
    # Store the lengths of the sequences
    n_i, n_j = len(seq_i), len(seq_j)
    if n_i > n_j:
        n_i, n_j = n_j, n_i
        seq_i, seq_j = seq_j, seq_i
    # Initialize scores and traceback matrices
    def scores_arrows_matrix(a, b, gap):
        mat_scores = [[0 for _ in range(n_j + 1)] for _ in range(n_i + 1)] # +1 for the initialization row/column
        mat_arrows = [[0 for _ in range(n_j + 1)] for _ in range(n_i + 1)] # +1 for the initialization row/column
        
        # Initialize edges with gaps
        for i in range(1, n_i + 1):
            mat_scores[i][0] = i * gap
            mat_arrows[i][0] = 1  # deletion (gap in seq_j)

        for j in range(1, n_j + 1):
            mat_scores[0][j] = j * gap
            mat_arrows[0][j] = -1  # insertion (gap in seq_i)
        
        for i in range(1, n_i + 1):
            for j in range(1, n_j + 1):
                diag = mat_scores[i-1][j-1] + subst_mat[str(a[i-1]), str(b[j-1])]
                left = mat_scores[i][j-1] + gap # left (insertion)
                up = mat_scores[i-1][j] + gap  # up (deletion)
                
                mat_scores[i][j] = max(diag, left, up)

                if diag > left and diag > up:
                    mat_arrows[i][j] = 0  # diagonal (match/mismatch)
                elif left > up:
                    mat_arrows[i][j] = -1  # left (insertion)
                else:
                    mat_arrows[i][j] = 1  # up (deletion)
        return mat_scores, mat_arrows
    scores, arrows = scores_arrows_matrix(seq_i, seq_j, gap) # seq_i = Rows / seq_j = Columns

    # Print scores matrix (for debugging)
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print(header)
    for aa, row in zip(f" {seq_i}", scores):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))
    
    # Print traceback matrix (for debugging))
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print("\n" + header)
    for aa, row in zip(f" {seq_i}", arrows):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

    # Print optimal score
    print(f"\nOptimal score: {scores[-1][-1]}") 

    # Prepare for traceback
    aln_i, aln_j = [], []
    i, j = len(arrows) - 1, len(arrows[0]) - 1 
    # Traceback
    while i > 0 or j > 0:  # fill in
        if arrows[i][j] == 0: # diagonal
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
        elif arrows[i][j] == -1:  # left
            j -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j])
        else:  # up (traceback[i][j] == 1)
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")
    # Print the alignment
    print("\nAlignment:")
    print("".join(aln_i[::-1]))
    print("".join(aln_j[::-1]))

nw('FAST', 'FAT', -2)