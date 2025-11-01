def needleman_wunsch_diagonal_loop(seq1, seq2, match=1, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize borders
    for i in range(1, m + 1):
        dp[i][0] = i * gap
    for j in range(1, n + 1):
        dp[0][j] = j * gap

    # Diagonal iteration
    for d in range(1, m + n):  # skip d=0 since that's just dp[0][0]
        for i in range(max(1, d - n + 1), min(m + 1, d + 1)):
            j = d - i
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            score = match if char1 == char2 else mismatch

            dp[i][j] = max(
                dp[i - 1][j - 1] + score,  # match/mismatch
                dp[i - 1][j] + gap,  # deletion
                dp[i][j - 1] + gap,  # insertion
            )

    return dp[m][n], dp


def needleman_wunsch_banded(seq1, seq2, bandwidth=3, match=1, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    INF = float("-inf")

    # Initialize DP table with -inf outside band
    dp = [[INF] * (n + 1) for _ in range(m + 1)]
    dp[0][0] = 0

    for i in range(1, m + 1):
        if i <= bandwidth:
            dp[i][0] = i * gap
    for j in range(1, n + 1):
        if j <= bandwidth:
            dp[0][j] = j * gap

    # Antidiagonal traversal with band constraint: |i - j| <= bandwidth
    for d in range(1, m + n):
        for i in range(max(1, d - n + 1), min(m + 1, d + 1)):
            j = d - i
            if abs(i - j) > bandwidth:
                continue

            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            score = match if char1 == char2 else mismatch

            best = INF
            if abs(i - 1 - j) <= bandwidth:
                best = max(best, dp[i - 1][j] + gap)
            if abs(i - j + 1) <= bandwidth:
                best = max(best, dp[i][j - 1] + gap)
            if abs(i - 1 - (j - 1)) <= bandwidth:
                best = max(best, dp[i - 1][j - 1] + score)

            dp[i][j] = best

    return dp[m][n] if abs(m - n) <= bandwidth else None, dp


def needleman_wunsch_banded_loops(
    seq1, seq2, bandwidth=3, match=1, mismatch=-1, gap=-1
):
    m, n = len(seq1), len(seq2)
    INF = float("-inf")

    # Initialize DP table with -inf
    dp = [[INF] * (n + 1) for _ in range(m + 1)]
    dp[0][0] = 0

    # Fill first row and column within the band
    for i in range(1, m + 1):
        if i <= bandwidth:
            dp[i][0] = i * gap
    for j in range(1, n + 1):
        if j <= bandwidth:
            dp[0][j] = j * gap

    for d in range(1, m + n):  # antidiagonals from i+j = 1 to m+n-1
        # Loop over i s.t. i + j = d and |i - j| <= bandwidth
        i_min = max(1, d - n + 1, (d - bandwidth + 1) // 2)
        i_max = min(m, d, (d + bandwidth) // 2)

        for i in range(i_min, i_max + 1):
            j = d - i
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            score = match if char1 == char2 else mismatch
            best = INF
            if i > 0 and dp[i - 1][j] != INF:
                best = max(best, dp[i - 1][j] + gap)
            if j > 0 and dp[i][j - 1] != INF:
                best = max(best, dp[i][j - 1] + gap)
            if i > 0 and j > 0 and dp[i - 1][j - 1] != INF:
                best = max(best, dp[i - 1][j - 1] + score)

            dp[i][j] = best

    return dp[m][n] if abs(m - n) <= bandwidth else None, dp


if __name__ == "__main__":
    seq1 = "A" * 3000
    seq2 = "A" * 3000
    # time each function
    import time

    start_time = time.time()
    score_norm, dp_matrix_norm = needleman_wunsch_diagonal_loop(seq1, seq2)
    end_time = time.time()
    print("Normal Needleman-Wunsch Time:", end_time - start_time)
    start_time = time.time()
    score_band, dp_matrix_band = needleman_wunsch_banded(seq1, seq2, bandwidth=2)
    end_time = time.time()
    print("Banded Needleman-Wunsch Time:", end_time - start_time)
    start_time = time.time()
    score_loop_band, dp_matrix_loop_band = needleman_wunsch_banded_loops(
        seq1, seq2, bandwidth=2
    )
    end_time = time.time()
    print("Loop Banded Needleman-Wunsch Time:", end_time - start_time)
    # print score norm
    # print score band
    # print score loop band
#     print("Score (Normal):", score_norm)
#     print("Score (Banded):", score_band)
#     print("Score (Loop Banded):", score_loop_band)
#     print("DP Matrix (Normal):")
#     for row in dp_matrix_norm:
#         print(row)
#     print("DP Matrix (Banded):")
#     for row in dp_matrix_band:
#         print(row)
#     print("DP Matrix (Loop Banded):")
#     for row in dp_matrix_loop_band:
#         print(row)
# # This code implements the Needleman-Wunsch algorithm with diagonal and banded optimizations.
