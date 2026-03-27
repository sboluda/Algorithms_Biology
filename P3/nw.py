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
    def scores_matrix(a, b, match, mismatch, gap):
        def scores_matrix_compute(a, b, mat, match, mismatch, gap):
            for i in range(len(a)):
                for j in range(len(b)):
                    if i == 0 and j == 0:
                        pass
                    elif i == 0:
                        mat[i][j] = mat[i][j-1] + gap
                    elif j == 0:
                        mat[i][j] = mat[i-1][j] + gap
                    else:
                        matching = 0
                        if a[i] == b[j]:
                            matching = mat[i-1][j-1] + match
                        else:
                            matching = mat[i-1][j-1] + mismatch
                        insertion = mat[i-1][j] + gap
                        deletion = mat[i][j-1] + gap
                        mat[i][j] = max(matching, insertion, deletion)
            return mat

        mat = [[0 for _ in range(len(b) + 1)] for _ in range(len(a) + 1)] # +1 for the initialization row/column
        return scores_matrix_compute("-" + a, "-" + b, mat, match, mismatch, gap) # Adding a dash for the initializations row/column
    scores = scores_matrix(seq_i, seq_j, match, mismatch, gap) # seq_i = Rows / seq_j = Columns

    # Print scores matrix (for debugging)
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print(header)
    for aa, row in zip(f" {seq_i}", scores):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

    arrow_scores = {"up": 1, "left": -1, "diag": 0}
    def traceback_matrix(a, b, scores):
        def traceback_matrix_compute(a, b, scores, mat):
            for i in range(len(a)):
                for j in range(len(b)):
                    if i == 0 and j == 0:
                        mat[i][j] == 0
                    elif i == 0:
                        mat[i][j] = -1
                    elif j == 0:
                        mat[i][j] = 1
                    else:
                        # Map varaibles names to values
                        if a[i] == b[j]:
                            diag = scores[i-1][j-1] + match
                        else:
                            diag = scores[i-1][j-1] + mismatch

                        up = scores[i-1][j] + gap
                        left = scores[i][j-1] + gap

                        possibilities = {
                            "diag": diag,
                            "up": up,
                            "left": left
                        }
                         # Fill the score matrix
                        scores[i][j] = max(diag, left, up)

                        # Fill the traceback matrix deterministically
                        if diag >= left and diag >= up:
                            mat[i][j] = 0  # diagonal
                        elif left >= up:
                            mat[i][j] = -1  # left
                        else:
                            mat[i][j] = 1  # up

            mat[-1][-1] = 0 # Set the last cell as 0 (where we will start)
            return mat

        mat = [[0 for _ in range(len(b) + 1)] for _ in range(len(a) + 1)] # +1 for the initialization row/column
        return traceback_matrix_compute("-" + a, "-" + b, scores, mat) # Adding a dash for the initializations row/column
    traceback = traceback_matrix(seq_i, seq_j, scores) # seq_i = Rows / seq_j = Columns

    # Print traceback matrix (for debugging))
    header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
    print("\n" + header)
    for aa, row in zip(f" {seq_i}", traceback):
        print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

    # Print optimal score
    print(f"Optimal score: {scores[-1][-1]}") 

    # Prepare for traceback
    aln_i, aln_j = [], []
    i, j = len(traceback) - 1, len(traceback[0]) - 1 
    # Traceback
    while i > 0 or j > 0:  # fill in
        if traceback[i][j] == 0: # diagonal
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
        elif traceback[i][j] == -1:  # left
            j -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j])
        else:  # up (traceback[i][j] == 1)
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")
    # Print the alignment
    print("Alignment:")
    print("".join(aln_i[::-1]))
    print("".join(aln_j[::-1]))

# seq_i = "FAST"
# seq_j = "FAT"
# match = 2 # Actual should be 1
# mismatch = -1 
# gap = -2 # Actual should be -2

nw(args.seq_i, args.seq_j, args.match, args.mismatch, args.gap)


