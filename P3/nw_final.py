import argparse

parser = argparse.ArgumentParser()
parser.add_argument("seq_i", type=str, help="First sequence")
parser.add_argument("seq_j", type=str, help="Second sequence")
parser.add_argument("match", type=int, help="Match score")
parser.add_argument("mismatch", type=int, help="Mismatch score")
parser.add_argument("gap", type=int, help="Gap score")

args = parser.parse_args()

def nw(seq_i, seq_j, match, mismatch, gap):
    # Store the lengths of the sequences
    n_i, n_j = len(seq_i), len(seq_j)

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

                if diag > left and diag > up:
                    mat_arrows[i][j] = 0  # diagonal (match/mismatch)
                elif left > up:
                    mat_arrows[i][j] = -1  # left (insertion)
                else:
                    mat_arrows[i][j] = 1  # up (deletion)
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
    print("\nAlignment:")
    print("".join(aln_i[::-1]))
    print("".join(aln_j[::-1]))

# seq_i = "FAST"
# seq_j = "FAT"
# match = 2 # Actual should be 1
# mismatch = -1 
# gap = -2 # Actual should be -2

nw(args.seq_i, args.seq_j, args.match, args.mismatch, args.gap)

#-------------------------------------------------------------------------#
### Exercise 2: Which parameter you have to decrease (match, mismatch or gap) to produce the following alignment? ###
# THE-BIGCAT
# THER----AT

# Reducing the Gap value to 0 allows us to get that alignment, given that the first - after "THE" implies a gap being better than a mismatch

### Exercise 3: Alternatively, you could have obtained the alignment in the same alignment above by letting the previous parameter untouched and increasing another parameter (match, mismatch or gap). Which one?

# In this case, if we said that by decreasing the gap penalty we favour it rather than the mismatch, the alternative is to worsen the mismatch. By just giving a mismatch value of -4 we can get that

### Exercise 4: Re-align the sequences THEBIGCAT and THERAT using a match score of +8, a mismatch score of −8, and a gap penalty of −4. Then reverse the order in which the sequences are provided to the program and recompute the alignment. What do you observe, and why does this happen?

# We observe...
# Optimal score: 20
# Alignment:
# THE----RAT
# THEBIGC-AT

# 1. The Optimal score is the same as if the order of the sequences was as before
# 2. The gaps are still the same. As we see gap_value > mismatch_value, hence we prefer gaps over mismatches
# 3. The align difference comes from the sequence positioning and how are the values located in the matrix. In this case, gaps appear earlier in sequence 1 (THERAT), while in the normal way, gaps appear first in the sequence 2 (THEBIGCAT)
# TL:DR, Same score, different path
