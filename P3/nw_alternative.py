import argparse

parser = argparse.ArgumentParser()
parser.add_argument("seq_i", type=str, help="First sequence")
parser.add_argument("seq_j", type=str, help="Second sequence")
parser.add_argument("match", type=int, help="Match score")
parser.add_argument("mismatch", type=int, help="Mismatch score")
parser.add_argument("gap", type=int, help="Gap score")

args = parser.parse_args()

seq_i = args.seq_i
seq_j = args.seq_j
match = args.match
mismatch = args.mismatch
gap = args.gap

# Store the lengths of the sequences
n_i, n_j = len(seq_i), len(seq_j)

# Initialize scores and traceback matrices
scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

# Initialize edges with gaps
for i in range(1, n_i + 1):
    scores[i][0] = i * gap
    traceback[i][0] = 1  # deletion (gap in seq_j)

for j in range(1, n_j + 1):
    scores[0][j] = j * gap
    traceback[0][j] = -1  # insertion (gap in seq_i)

# Fill the matrices
for i in range(1, n_i + 1):
    for j in range(1, n_j + 1):
        # Calculate scores for all possibilities
        match_score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch
        diag = scores[i - 1][j - 1] + match_score  # match/mismatch
        left = scores[i][j - 1] + gap  # insertion (gap in seq_i)
        up = scores[i - 1][j] + gap  # deletion (gap in seq_j)

        # Choose the best scores
        scores[i][j] = max(diag, left, up)

        # Set traceback pointer
        if diag > left and diag > up:
            traceback[i][j] = 0  # diagonal (match/mismatch)
        elif left > up:
            traceback[i][j] = -1  # left (insertion)
        else:
            traceback[i][j] = 1  # up (deletion)

# Print scores matrix
header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
print(header)
for aa, row in zip(f" {seq_i}", scores):
    print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

# Print traceback matrix
header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
print("\n" + header)
for aa, row in zip(f" {seq_i}", traceback):
    print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

# Print optimal score
print("\nOptimal score:", scores[-1][-1], "\n")

# Traceback
aln_i, aln_j = [], []
i, j = n_i, n_j

while i > 0 or j > 0:
    if traceback[i][j] == 0:  # diagonal
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
print("".join(reversed(aln_i)))
print("".join(reversed(aln_j)))