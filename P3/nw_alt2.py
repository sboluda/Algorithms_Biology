import argparse

parser = argparse.ArgumentParser()
parser.add_argument("seq_i", type=str, help="First sequence")
parser.add_argument("seq_j", type=str, help="Second sequence")
parser.add_argument("match", type=int, help="Match score")
parser.add_argument("mismatch", type=int, help="Mismatch score")
parser.add_argument("gap", type=int, help="Gap score")

args = parser.parse_args()

def nw(seq_i, seq_j, match, mismatch, gap): # We change the order of the parameters!!
    # Store the lengths of the sequences
    n_i, n_j = len(seq_i), len(seq_j)
    if n_i < n_j:
        n_i, n_j = n_j, n_i
        seq_i, seq_j = seq_j, seq_i

    # Initialize scores and traceback matrices
    def scores_arrows_matrix(a, b, match, mismatch, gap):
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
                diag = mat_scores[i-1][j-1] + (match if a[i-1]==b[j-1] else mismatch)
                left = mat_scores[i][j-1] + gap # left (insertion)
                up = mat_scores[i-1][j] + gap  # up (deletion)
                
                mat_scores[i][j] = max(diag, left, up)

                if up > left and up > diag: # We check up first!!
                    mat_arrows[i][j] = 1  # up (deletion)
                elif left > up:
                    mat_arrows[i][j] = -1  # left (insertion)
                else:                      # We check diagonal last (tie case)!!
                    mat_arrows[i][j] = 0  # diagonal (match/mismatch)
        return mat_scores, mat_arrows
    scores, arrows = scores_arrows_matrix(seq_i, seq_j, match, mismatch, gap) # seq_i = Rows / seq_j = Columns

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
    print("\nAlignment:") # We change the order of the prints here!!
    print("".join(aln_i[::-1]))
    print("".join(aln_j[::-1])) 

nw(args.seq_i, args.seq_j, args.match, args.mismatch, args.gap)

# The only change has been in the create of the traceback/arrows matrix. Which considered that, in case of a tie, to output the +1 (up)

# Now we just need to change it so that in all those ties Diagonal (Match/Mismatch) is the default case, that way we prioritize that sceneario