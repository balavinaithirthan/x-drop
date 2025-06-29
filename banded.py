def needleman_wunsch_banded_loops(
    seq1, seq2, bandwidth=3, match=1, mismatch=-1, gap=-1
):
    m, n = len(seq1), len(seq2)
    INF = float("-inf")

    # Initialize DP table with -inf
    dp = [[INF] * (n + 1) for _ in range(m + 1)]
    dp[0][0] = 0

    # Fill first row and column within the band
    # for i in range(1, m + 1):
    #     if i <= bandwidth:
    #         dp[i][0] = i * gap
    # for j in range(1, n + 1):
    #     if j <= bandwidth:
    #         dp[0][j] = j * gap

    for d in range(0, m + n + 1):  # antidiagonals from i+j = 1 to m+n-1
        # Loop over i s.t. i + j = d and |i - j| <= bandwidth
        i_min = max(0, d - n, (d - bandwidth + 1) // 2)
        i_max = min(m, d, (d + bandwidth) // 2)
        # i_min = max(0, d - n)
        # i_max = min(m, d)

        for i in range(i_min, i_max + 1):
            j = d - i
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            score = match if char1 == char2 else mismatch

            # best = INF
            # if i > 0 and dp[i - 1][j] != INF:
            #     best = max(best, dp[i - 1][j] + gap)
            # if j > 0 and dp[i][j - 1] != INF:
            #     best = max(best, dp[i][j - 1] + gap)
            # if i > 0 and j > 0 and dp[i - 1][j - 1] != INF:
            #     best = max(best, dp[i - 1][j - 1] + score)

            dp[i][j] = 10

    return dp[m][n] if abs(m - n) <= bandwidth else None, dp


def generated_needleman_wunsch_banded(
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
        #    for (int ad = 0; ad < ((m + n) + 1); ++ad) {
        #     for (int i = (max((((ad / 2) + (3 / 2))), max(0, (ad - m)))); i < min((min(n, ad) + 1), (((ad / 2) - (3 / 2)))); ++i) {
        #         j = (ad - i);
        #     }
        # }
        #
        # i_min = max(1, d - n + 1, (d - bandwidth + 1) // 2)
        # i_max = min(m, d, (d + bandwidth) // 2)
        for i in range(
            max((d // 2) - (3 // 2), 1, d - n),
            min(min(m, d) + 1, (d // 2) + (3 // 2)),
        ):
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


seq1 = "A" * 10
seq2 = "A" * 10
score_loop_band, dp_matrix_loop_band = needleman_wunsch_banded_loops(
    seq1, seq2, bandwidth=100
)
# score_loop_band_gen, dp_matrix_loop_band_gen = generated_needleman_wunsch_banded(
# seq1, seq2, bandwidth=15
# )

for row in dp_matrix_loop_band:
    print(row)

    # print("Score (generated banded):", score_loop_band_gen)
    # for row in dp_matrix_loop_band_gen:
    print(row)
