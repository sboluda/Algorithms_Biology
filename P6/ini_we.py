from Bio.Align import substitution_matrices


def we(seq_i, seq_j, gap, n_alignments=5, threshold=0):
    # sourcery skip: low-code-quality
    """
    Waterman-Eggert implementation to find multiple optimal local alignments.

    Args:
        seq_i: First sequence
        seq_j: Second sequence
        gap: Gap penalty
        n_alignments: Maximum number of alignments to return
        threshold: Minimum score to consider an alignment

    Returns:
        List of tuples (score, (aln_i, aln_j, start_i, end_i, start_j, end_j))

    >>> alignments = we("FAT", "THEFASTCAT", -4, n_alignments=3, threshold=0)
    >>> len(alignments)
    3
    >>> alignments[0][0]
    11.0
    >>> alignments[0][1]
    ('FAT', 'FAS', 0, 3, 3, 6)
    >>> alignments[0][1][:2]
    ('FAT', 'FAS')
    >>> alignments[1][0]
    9.0
    >>> alignments[2][0]
    6.0
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    alignments = []

    # Matrix to mark cells used in previous alignments

    used = [[]]  # FILL IN: set to False for all cells

    # Loop to find multiple alignments 
    for _ in range(n_alignments): 

        # Build the score matrix
        scores = [[0 for _ in range(len(seq_j) + 1)] for _ in range(len(seq_i) + 1)]
        traceback = [[2 for _ in range(len(seq_j) + 1)] for _ in range(len(seq_i) + 1)]

        best_score = 0
        best_score_cell = (0, 0)

        for i in range(1, len(seq_i) + 1):
            for j in range(1, len(seq_j) + 1):
                # If the cell was already used, set score to 0 and the traceback to 2

                # FILL IN
                
                    continue

                # Calculate scores for all possibilities
                match_score = subst_mat[seq_i[i - 1], seq_j[j - 1]]
                diag = scores[i - 1][j - 1] + match_score  # match/mismatch
                left = scores[i][j - 1] + gap  # insertion (gap in seq_i)
                up = scores[i - 1][j] + gap  # deletion (gap in seq_j)

                # Choose the best score
                scores[i][j] = max(0, diag, left, up)

                # Set traceback pointer
                if diag > left and diag > up and diag > 0:
                    traceback[i][j] = 0  # diagonal (match/mismatch)
                elif left > up and left > 0:
                    traceback[i][j] = -1  # left (insertion)
                elif up > 0:
                    traceback[i][j] = 1  # up (deletion)
                else:
                    traceback[i][j] = 2

                if scores[i][j] > best_score:
                    best_score = scores[i][j]
                    best_score_cell = (i, j)

        # Stop if best score below threshold
        if best_score < threshold:
            break

        # Traceback and mark used cells
        i, j = best_score_cell
        end_i, end_j = i, j
        aln_i = []
        aln_j = []
        path = []

        while traceback[i][j] != 2:
            path.append((i, j))
            if traceback[i][j] == 0:  # Match/mismatch
                i -= 1
                j -= 1
                aln_i.append(seq_i[i])
                aln_j.append(seq_j[j])
            elif traceback[i][j] == -1:  # Deletion
                j -= 1
                aln_i.append("-")
                aln_j.append(seq_j[j])
            else:  # Insertion
                i -= 1
                aln_i.append(seq_i[i])
                aln_j.append("-")

        start_i, start_j = i, j

        # Mark all cells in the path as used

        # FILL IN

        seq_i_aln = "".join(reversed(aln_i))
        seq_j_aln = "".join(reversed(aln_j))

        alignments.append(
            (best_score, (seq_i_aln, seq_j_aln, start_i, end_i, start_j, end_j))
        )

    return alignments


if __name__ == "__main__":
    print("=" * 60)
    print("Waterman-Eggert (multiple alignments)")
    print("=" * 60)
    alignments = we("FAT", "THEFASTCAT", -4, n_alignments=3, threshold=0)

    for idx, (score, (aln_i, aln_j, start_i, end_i, start_j, end_j)) in enumerate(
        alignments, 1
    ):
        print(f"Alignment {idx}:")
        print(f"  Score: {score}")
        print(f"  seq_i [{start_i}:{end_i}]:\t{aln_i}")
        print(f"  seq_j [{start_j}:{end_j}]:\t{aln_j}")
        print()

    print("=" * 60)
    print("Example with higher similarity")
    print("=" * 60)
    alignments = we("THECATISFAT", "AFASTCAT", -2, n_alignments=5, threshold=6)

    for idx, (score, (aln_i, aln_j, start_i, end_i, start_j, end_j)) in enumerate(
        alignments, 1
    ):
        print(f"Alignment {idx}:")
        print(f"  Score: {score}")
        print(f"  seq_i [{start_i}:{end_i}]:\t{aln_i}")
        print(f"  seq_j [{start_j}:{end_j}]:\t{aln_j}")
        print()
