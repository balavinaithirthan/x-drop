from math import inf

def x_drop(seq1, seq2, x_drop=100, match=0, mismatch=1, gap=1):
    """
    Simplified x-drop alignment scoring function using antidiagonals.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        match (int): Score for match
        mismatch (int): Score for mismatch
        gap (int): Score for gap

    Returns:
        tuple: (final score, DP matrix)
    """
    m, n = len(seq1), len(seq2)
    dp = [[inf for _ in range(n + 1)] for _ in range(m + 1)]

    var_lo = inf
    var_hi = -inf

    for d in range(m + n + 1):  # iterate over antidiagonals
        i_start = max(0, d - n)
        i_start = max(i_start, (d - var_hi) // 2 )

        i_end = min(m, d)
        i_end = min(i_end, (d + var_lo) // 2)

        for i in range(i_start, i_end + 1):
            j = d - i

            # Initialize first row/col with gap penalties
            if i == 0:
                dp[i][j] = j * gap
            elif j == 0:
                dp[i][j] = i * gap
            else:
                diag = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                up = dp[i-1][j] + gap
                left = dp[i][j-1] + gap
                dp[i][j] = min(diag, up, left)
        
            if dp[i][j] > x_drop:
                dp[i][j] = inf

        # Update pruning variables (not crucial here but retained from original)
        # next_i = None
        # for i in range(i_end, i_start - 1, -1):
        #     j = d - i
        #     if dp[i][j] != inf:
        #         next_i = i
        #         break
        # if next_i is not None and next_i != i_end:
        #     j = d - next_i
        #     var_lo = min(j - next_i, var_lo)

        # next_i = None
        # for i in range(i_start, i_end + 1):
        #     j = d - i
        #     if dp[i][j] != inf:
        #         next_i = i
        #         print("next_i is:", next_i)
        #         break
        # if next_i is not None and next_i != i_start:
        #     j = d - next_i
        #     var_hi = max(j - next_i, var_hi)
        
        # if var_lo > var_hi:
        #     return dp[m][n], dp
        
        # print(var_lo, var_hi)

    return min(dp[m][n], x_drop), dp


def print_dp_matrix(dp, seq1, seq2):
    """Nicely prints the DP matrix with sequence labels."""
    print("     ", end="")
    for c in " " + seq2:
        print(f"{c:>4}", end="")
    print()
    for i, row in enumerate(dp):
        label = seq1[i-1] if i > 0 else " "
        print(f"{label:>3} ", end="")
        for val in row:
            print(f"{val:>4}", end="")
        print()
    print()


# ------------------- TEST CASES -------------------

def run_tests():
    test_cases = [
        ("GATT", "GCT"),
        ("AAGGCG", "ACGGCG"),
        ("AAA", "AA"),
        ("ACGT", "ACGT"),
        ("GCGT", "GAT"),
    ]

    for idx, (seq1, seq2) in enumerate(test_cases, start=1):
        print(f"=== Test {idx}: seq1='{seq1}', seq2='{seq2}' ===")
        score, dp = x_drop(seq1, seq2)
        print(f"Final alignment score: {score}")
        print_dp_matrix(dp, seq1, seq2)


if __name__ == "__main__":
    run_tests()
