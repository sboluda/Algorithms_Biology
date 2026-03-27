def nw_score_only(seq_i, seq_j, match, mismatch, gap):
    """Computes exactly the same scores as NW, but discards everything except the last
    row, making traceback impossible while allowing linear memory usage.
    Returns only the last row of the NW DP matrix.
    >>> nw_score_only ("FAT", "FAST", 2, -1, -1)
    [-3, 0, 3, 3, 5]
    """
    prev = [j * gap for j in range(len(seq_j) + 1)]

    for i in range(1, len(seq_i) + 1):
        curr = [i * gap]
        for j in range(1, len(seq_j) + 1):
            match_score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch
            diag = prev[j - 1] + match_score
            left = curr[j - 1] + gap
            up = prev[j] + gap
            curr.append(max(diag, left, up))
        prev = curr

    return prev


def nw_basic(seq_i, seq_j, match, mismatch, gap):
    """Classic Needleman-Wunsch (used as base case).
    >>> nw_score_only ("FAT", "FAST", 2, -1, -1)
    [-3, 0, 3, 3, 5]
    """
    n_i, n_j = len(seq_i), len(seq_j)

    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    for i in range(1, n_i + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1  # up (gap in seq_j)

    for j in range(1, n_j + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1  # left (gap in seq_i)

    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            match_score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch
            diag = scores[i - 1][j - 1] + match_score
            left = scores[i][j - 1] + gap
            up = scores[i - 1][j] + gap

            scores[i][j] = max(diag, left, up)

            if diag >= left and diag >= up:
                traceback[i][j] = 0  # diagonal (match/mismatch)
            elif left >= up:
                traceback[i][j] = -1  # left (gap in seq_i)
            else:
                traceback[i][j] = 1  # up (gap in seq_j)

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
        else:  # up
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")

    return "".join(reversed(aln_i)), "".join(reversed(aln_j))


def myers_miller(seq_i, seq_j, match, mismatch, gap, cutoff=2):
    """
    Global alignment using Myers-Miller (linear-space Needleman-Wunsch).

    >>> myers_miller("FAT", "FAST", 2, -1, -1)
    ('FA-T', 'FAST')
    >>> myers_miller("THEFASTCAT", "THERAT", 1, -8, -4)
    ('THE-FASTCAT', 'THER-A-T---')
    >>> myers_miller("THERAT", "THEFASTCAT", 1, -8, -4)
    ('THE-RA-T---', 'THEF-ASTCAT')
    """
    # Base case
    if len(seq_i) <= cutoff or len(seq_j) <= cutoff:
        return nw_basic(seq_i, seq_j, match, mismatch, gap)

    mid = len(seq_i) // 2

    score_forward = nw_score_only(seq_i[:mid], seq_j, match, mismatch, gap)
    score_backward = nw_score_only(seq_i[mid:][::-1], seq_j[::-1], match, mismatch, gap)

    # Find the column of seq_j through which the optimal alignment passes at row mid
    best_j = 0
    best_score = float("-inf")
    for j in range(len(seq_j) + 1):
        score = score_forward[j] + score_backward[len(seq_j) - j]
        if score > best_score:
            best_score = score
            best_j = j

    aln_i_forward, aln_j_forward = myers_miller(
        seq_i[:mid], seq_j[:best_j], match, mismatch, gap, cutoff
    )
    aln_i_backward, aln_j_backward = myers_miller(
        seq_i[mid:], seq_j[best_j:], match, mismatch, gap, cutoff
    )

    return aln_i_forward + aln_i_backward, aln_j_forward + aln_j_backward


if __name__ == "__main__":
    aln_i, aln_j = myers_miller("THEFASTCAT", "THERAT", 1, 0, -1)
    print(aln_i)
    print(aln_j)
