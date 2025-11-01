def edit_distance(seq1, seq2, Z=5, match=0, mismatch=-1, gap=-1):
    """
    Standard dynamic programming edit distance (Levenshtein distance).

    Args:
        seq1, seq2 (str): Input sequences.
        Z (int): Unused placeholder (for consistency).
        match, mismatch, gap (int): Scoring parameters (not used directly here).

    Returns:
        D (list[list[int]]): Full DP table.
        dist (int): Minimum edit distance.
    """
    m, n = len(seq1), len(seq2)
    D = [[0] * (n + 1) for _ in range(m + 1)]
    cells_explored = 0

    # Initialization
    for i in range(m + 1):
        D[i][0] = i
        cells_explored += 1
    for j in range(1, n + 1):  # Skip [0][0] already counted
        D[0][j] = j
        cells_explored += 1

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            D[i][j] = min(
                D[i - 1][j] + 1,      # deletion
                D[i][j - 1] + 1,      # insertion
                D[i - 1][j - 1] + cost  # substitution
            )
            cells_explored += 1

    def debug_print(D):
        for i in range(len(D)):
            for j in range(len(D[i])):
                print(D[i][j], end=" ")
            print()
        print()
    # debug_print(D)


    return D, D[m][n], cells_explored
