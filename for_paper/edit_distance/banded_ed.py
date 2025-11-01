from math import inf

def banded_edit_distance(seq1, seq2, bandwidth=30):
    """
    Banded edit distance DP: computes within a diagonal band of width `bandwidth`.

    Args:
        seq1, seq2 (str): Input sequences.
        bandwidth (int): Half-width of the diagonal band.
        Z (int): Placeholder for consistency.
        match, mismatch, gap (int): Scoring parameters (not used directly here).

    Returns:
        D (list[list[int]]): Banded DP table (incomplete outside the band).
        dist (int): Minimum edit distance if within the band; otherwise approximate.
        cells_explored (int): Number of DP cells actually computed.
    """
    m, n = len(seq1), len(seq2)
    D = [[inf] * (n + 1) for _ in range(m + 1)]
    cells_explored = 0

    # Initialization
    D[0][0] = 0
    cells_explored += 1
    for i in range(1, m + 1):
        if abs(i - 0) <= bandwidth:
            D[i][0] = i
            cells_explored += 1
    for j in range(1, n + 1):
        if abs(0 - j) <= bandwidth:
            D[0][j] = j
            cells_explored += 1

    # Fill only inside the band
    for i in range(1, m + 1):
        j_min = max(1, i - bandwidth)
        j_max = min(n, i + bandwidth)
        for j in range(j_min, j_max + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            D[i][j] = min(
                D[i - 1][j] + 1,      # deletion
                D[i][j - 1] + 1,      # insertion
                D[i - 1][j - 1] + cost  # substitution
            )
            cells_explored += 1

    return D, D[m][n], cells_explored


# --------------------------- TESTING SCRIPT ---------------------------

def print_matrix(D, seq1, seq2, bandwidth):
    """Prints DP matrix, showing INF as '∞' for clarity."""
    print(f"\nDP Matrix (bandwidth={bandwidth}):")
    print("     ", end="")
    for c in " " + seq2:
        print(f"{c:>4}", end="")
    print()
    for i, row in enumerate(D):
        label = seq1[i - 1] if i > 0 else " "
        print(f"{label:>3} ", end="")
        for val in row:
            print(f"{'∞' if val > 1000 else val:>4}", end="")
        print()
    print()


def run_tests():
    tests = [
        ("A" * 5 + "B" * 5, "B" * 5 + "A" * 5, 5),
        ("GATTACA", "GCATGCU", 2),
        ("ACGT", "ACGT", 1),
        ("AAAA", "TTTT", 2),
        ("ABCDEF", "AZCED", 3),
        ("SUNDAY", "SATURDAY", 3),
        ("HELLO", "HE", 1)
    ]

    for idx, (s1, s2, bw) in enumerate(tests, start=1):
        print(f"=== Test {idx}: seq1='{s1}', seq2='{s2}', bandwidth={bw} ===")
        D, dist, explored = banded_edit_distance(s1, s2, bandwidth=bw)
        print(f"Edit distance: {dist}")
        print(f"Cells explored: {explored} / {(len(s1)+1)*(len(s2)+1)} total")
        print_matrix(D, s1, s2, bw)


if __name__ == "__main__":
    run_tests()
