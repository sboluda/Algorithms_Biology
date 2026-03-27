def gotoh(seq_i, seq_j, match, mismatch, gap_open, gap_extend):
    """
    Needleman-Wunsch global alignment with affine gap penalties (Gotoh Algorithm).

    seq_i is on the Y axis (rows, index i), seq_j is on the X axis (columns, index j).

    M[i][j]  = best score ending with seq_i[i-1] aligned to seq_j[j-1]
    Ix[i][j] = best score ending with a gap in seq_j (move up/vertically, consuming seq_i[i-1])
    Iy[i][j] = best score ending with a gap in seq_i (move left/horizontally, consuming seq_j[j-1])

    Traceback convention (consistent with NW):
        0  = diagonal (match/mismatch)
        1  = up       (gap in seq_j, Ix)
        -1 = left     (gap in seq_i, Iy)

    >>> gotoh("FAT", "FAST", 2, -1, -2, -1)
    ('FA-T', 'FAST')
    >>> gotoh("THEFASTCAT", "THERAT", 1, -3, -3, -1)
    ('THEFASTCAT', 'THE----RAT')
    """
    n, m = len(seq_i), len(seq_j)
    NEG_INF = float("-inf")

    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # tb_M:  0=from M, 1=from Ix, -1=from Iy  (all diagonal moves)
    # tb_Ix: 0=opened from M,      1=extended Ix  (all upward moves)
    # tb_Iy: 0=opened from M,     -1=extended Iy  (all leftward moves)
    tb_M = [[0] * (m + 1) for _ in range(n + 1)]
    tb_Ix = [[0] * (m + 1) for _ in range(n + 1)]
    tb_Iy = [[0] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0

    # First column: only Ix gaps are valid (moving up, gap in seq_j)
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + i * gap_extend
        tb_Ix[i][0] = 1  # extend Ix (up)

    # First row: only Iy gaps are valid (moving left, gap in seq_i)
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + j * gap_extend
        tb_Iy[0][j] = -1  # extend Iy (left)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch

            # M: best alignment ending in a match/mismatch (diagonal move)
            m_val = M[i - 1][j - 1] + score
            ix_val = Ix[i - 1][j - 1] + score
            iy_val = Iy[i - 1][j - 1] + score

            if m_val >= ix_val and m_val >= iy_val:
                M[i][j], tb_M[i][j] = m_val, 0
            elif ix_val >= iy_val:
                M[i][j], tb_M[i][j] = ix_val, 1
            else:
                M[i][j], tb_M[i][j] = iy_val, -1

            # Ix: gap in seq_j, move up (consume seq_i[i-1])
            open_ix = M[i - 1][j] + gap_open + gap_extend
            ext_ix = Ix[i - 1][j] + gap_extend

            if open_ix >= ext_ix:
                Ix[i][j], tb_Ix[i][j] = open_ix, 0
            else:
                Ix[i][j], tb_Ix[i][j] = ext_ix, 1

            # Iy: gap in seq_i, move left (consume seq_j[j-1])
            open_iy = M[i][j - 1] + gap_open + gap_extend
            ext_iy = Iy[i][j - 1] + gap_extend

            if open_iy >= ext_iy:
                Iy[i][j], tb_Iy[i][j] = open_iy, 0
            else:
                Iy[i][j], tb_Iy[i][j] = ext_iy, -1

    # Traceback: start from the best-scoring matrix at (n, m)
    aln_i, aln_j = [], []
    i, j = n, m

    if M[n][m] >= Ix[n][m] and M[n][m] >= Iy[n][m]:
        matrix = 0
    elif Ix[n][m] >= Iy[n][m]:
        matrix = 1
    else:
        matrix = 2

    while i > 0 or j > 0:
        if matrix == 0:  # M: aligned pair, diagonal move
            aln_i.append(seq_i[i - 1])
            aln_j.append(seq_j[j - 1])
            prev = tb_M[i][j]
            i -= 1
            j -= 1
            matrix = 0 if prev == 0 else (1 if prev == 1 else 2)
        elif matrix == 1:  # Ix: gap in seq_j, move up
            aln_i.append(seq_i[i - 1])
            aln_j.append("-")
            stays = tb_Ix[i][j] == 1
            i -= 1
            matrix = 1 if stays else 0
        else:  # Iy: gap in seq_i, move left
            aln_i.append("-")
            aln_j.append(seq_j[j - 1])
            stays = tb_Iy[i][j] == -1
            j -= 1
            matrix = 2 if stays else 0

    aln_i.reverse()
    aln_j.reverse()
    return "".join(aln_i), "".join(aln_j)


if __name__ == "__main__":
    print("Test 1:")
    aln = gotoh("FAT", "FAST", 2, -1, -2, -1)
    print(aln[0])
    print(aln[1])
    print("Expected: FA-T / FAST\n")

    print("Test 2:")
    aln = gotoh("THEFASTCAT", "THERAT", 1, -3, -3, -1)
    print(aln[0])
    print(aln[1])
    print("Expected: THEFASTCAT / THE----RAT")
