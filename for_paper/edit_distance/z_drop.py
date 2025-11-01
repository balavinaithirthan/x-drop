import math

def z_drop_edit_distance(seq1, seq2, Z=5, match=0, mismatch=1, gap=1):
    """
    Compute edit-distance-like alignment with Z-drop pruning along diagonals.
    seq1, seq2: strings
    Z: drop threshold; pruning occurs if (best - score) > Z + |Δgap|
    """
    m, n = len(seq1), len(seq2)
    dp = [[-math.inf] * (n + 1) for _ in range(m + 1)]

    # initialization
    dp[0][0] = 0
    for i in range(1, m + 1):
        dp[i][0] = i * gap
    for j in range(1, n + 1):
        dp[0][j] = j * gap

    best_score = 0
    best_i = best_j = 0
    diag_lo, diag_hi = 0, 0  # band boundaries
    cells_explored = 0

    # traverse anti-diagonals
    for d in range(1, m + n + 1):
        i_start = max(1, d - n)
        i_end = min(m, d)
        i_start = max(i_start, (d + diag_lo) // 2)
        i_end = min(i_end, (d - diag_hi) // 2)


        for i in range(i_start, i_end + 1):
            j = d - i
            if j < 1 or j > n:
                continue

            score_diag = dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            score_up   = dp[i - 1][j] + gap
            score_left = dp[i][j - 1] + gap
            score = max(score_diag, score_up, score_left)
            dp[i][j] = score
            cells_explored += 1

            if score > best_score:
                best_score = score
                best_i, best_j = i, j

            # z-drop diagonal pruning
            drop = best_score - score
            gap_penalty = abs((i - best_i) - (j - best_j))
            if drop > Z + gap_penalty:
                dp[i][j] = -math.inf

        # update band boundaries (diagonal range)
        active = [(d - i, i) for i in range(i_start, i_end + 1)
                  if 0 <= d - i < n + 1 and dp[i][d - i] != -math.inf]
        if active:
            diag_lo = min(diag_lo, active[0][0] - active[0][1])
            diag_hi = max(diag_hi, active[-1][0] - active[-1][1])
        
        # Early termination: if band has collapsed, stop computation
        if diag_lo >= diag_hi:
            return dp, dp[i][j], cells_explored, best_score, best_i, best_j, diag_lo, diag_hi

    return dp, dp[m][n], cells_explored, best_score, best_i, best_j, diag_lo, diag_hi


# ===================================================
# Test Harness
# ===================================================

def print_dp(seq1, seq2, dp):
    print(f"   {'  '.join('-' + c for c in ' ' + seq2)}")
    for i, row in enumerate(dp):
        prefix = '-' if i == 0 else seq1[i - 1]
        print(prefix, end="  ")
        print("  ".join("{:>4}".format("-inf" if x == -math.inf else int(x)) for x in row))
    print()


def run_test(seq1, seq2, Z):
    dp, final_score, cells_explored, best_score, bi, bj, lo, hi = z_drop_edit_distance(seq1, seq2, Z)
    print(f"=== Test: seq1='{seq1}' vs seq2='{seq2}' (Z={Z}) ===")
    print(f"Best score={best_score} at ({bi},{bj}); diag_lo={lo}, diag_hi={hi}")
    print_dp(seq1, seq2, dp)
    print(f"Final alignment score: {final_score}\n{'='*60}\n")


if __name__ == "__main__":
    # 1. Identical
    run_test("GATTACA", "GATTACA", Z=5)

    # 2. One mismatch
    run_test("GATTACA", "GACTACA", Z=5)

    # 3. Partial overlap with gap
    run_test("GATTACA", "GAATGCA", Z=4)

    # 4. Dissimilar
    run_test("AAAA", "TTTT", Z=5)

    # 5. Heavy pruning (small Z)
    run_test("GATTACA", "GCATGCU", Z=2)

    # 6. Minimal pruning (large Z)
    run_test("GATTACA", "GCATGCU", Z=10)
