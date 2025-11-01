def wfa_edit_distance(seq1, seq2, Z=5, match=0, mismatch=-1, gap=-1):
    """
    Compute edit distance using the Wavefront Alignment (WFA) algorithm.

    Args:
        seq1, seq2 (str): Input sequences.
        Z (int): Drop threshold (placeholder for z-drop pruning).
        match, mismatch, gap (int): Scoring parameters (currently unused).

    Returns:
        (W, e): Wavefront matrix and final edit distance.
    """

    P, T = seq1, seq2
    m, n = len(P), len(T)
    max_e = max(m, n)
    num_diags = 2 * max_e + 1
    offset = m
    diag_end = offset + (n - m)
    end_pt = n

    # Initialize wavefronts
    W = [[-1] * num_diags for _ in range(max_e + 1)]
    diag_lo = [0] * (max_e + 1)
    diag_hi = [0] * (max_e + 1)
    cells_explored = 0

    # Initial wavefront (e=0)
    W[0][offset] = 0
    cells_explored += 1
    v, h = 0, 0
    while v < m and h < n and P[v] == T[h]:
        v += 1
        h += 1
        W[0][offset] += 1
    diag_lo[0] = diag_hi[0] = 0

    # Perfect prefix match
    if W[0][diag_end] == end_pt:
        return W, 0, cells_explored

    # Wavefront expansion
    for e in range(1, max_e + 1):
        lo = diag_lo[e - 1] - 1
        hi = diag_hi[e - 1] + 1
        diag_lo[e] = lo
        diag_hi[e] = hi

        for k in range(lo, hi + 1):
            idx = offset + k
            if not (0 <= idx < num_diags):
                continue

            # Get predecessors safely
            from_del = W[e - 1][idx + 1] if 0 <= idx + 1 < num_diags else -1
            from_ins = W[e - 1][idx] + 1 if 0 <= idx < num_diags and W[e - 1][idx] >= 0 else -1
            from_sub = W[e - 1][idx - 1] + 1 if 0 <= idx - 1 < num_diags and W[e - 1][idx - 1] >= 0 else -1

            W[e][idx] = max(from_del, from_ins, from_sub)
            cells_explored += 1

            # Extend matches (guarded)
            i = W[e][idx]
            j = i - k
            while 0 <= i < m and 0 <= j < n and P[i] == T[j]:
                i += 1
                j += 1
                W[e][idx] += 1

            # Termination check
            if 0 <= diag_end < num_diags and W[e][diag_end] == end_pt:
                return W, e, cells_explored

    return W, max_e, cells_explored


if __name__ == "__main__":
    test_cases = [
    # Exact matches
    ("", "", 0),
    ("A", "A", 0),
    ("ACGT", "ACGT", 0),

    # Insertions / deletions
    ("A", "", 1),
    ("", "A", 1),
    ("ACGT", "ACG", 1),
    ("ACGT", "ACGTT", 1),

    # Substitutions
    ("A", "T", 1),
    ("ACGT", "AGGT", 1),
    ("ACGT", "TCGA", 2),

    # Complex examples
    ("GATTACA", "GCATGCU", 4),
    ("ACGTACGT", "", 8),
    ("", "ACGTACGT", 8),

    # Repeated patterns
    ("A" * 50, "A" * 50, 0),
    ("A" * 50, "T" * 50, 50),
    ("ACGT" * 20, "AGCT" * 20, 20),
    ]

    for seq1, seq2, expected in test_cases:
        _, e, cells = wfa_edit_distance(seq1, seq2)
        print(f"{seq1!r} vs {seq2!r} -> edit distance = {e} (expected {expected}), cells = {cells}")
