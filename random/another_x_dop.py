from math import inf


def create_banded_mask(m, n, bandwidth):

    matrix = [[-inf for _ in range(n)] for _ in range(m)]

    for i in range(m):
        for j in range(n):
            if abs(i - j) <= bandwidth:
                matrix[i][j] = 0
    return matrix


def print_matrix(matrix):
    for row in matrix:
        print("\t".join(f"{val:5}" if val != float("-inf") else " -inf" for val in row))


def x_drop(m, n, dp, bandwidth=2, match=1, mismatch=-1, gap=-1):
    # Fill first row and column within the band
    var_lo = inf
    var_hi = -inf
    for d in range(0, m + n + 1):  # antidiagonals from i+j = 1 to m+n-1
        # Loop over i s.t. i + j = d and |i - j| <= bandwidth
        # i_min = max(1, d - n + 1, (d - var_hi + 1) // 2)
        # i_max = min(m, d, (d + var_lo) // 2)
        i_start = max(0, d - n)
        i_start = max(i_start, (d + var_lo) // 2)
        i_end = min(m, d)
        i_end = min(i_end, (d - var_lo) // 2)
        for i in range(i_start, i_end + 1):
            j = d - i
            if d == 5:
                if abs(j - i) > 2:
                    dp[i][j] = -inf
                else:
                    dp[i][j] = 5
            else:
                if j == 0:
                    dp[i][0] = i * gap
                elif i == 0:
                    dp[0][j] = j * gap
                else:
                    dp[i][j] = 5

        # calculate i -> figure out next diagonal from that
        next_i = 0
        for i in range(i_end, i_start - 1, -1):
            j = d - i
            if dp[i][j] != -inf:
                next_i = i
                break
        # how to remove this guard
        if next_i != i_end:
            j = d - next_i
            var_lo = min(j - next_i, var_lo)  # prune axis

        next_i = 0
        for i in range(i_start, i_end):
            j = d - i
            if dp[i][j] != -inf:
                next_i = i
                break
        if next_i != i_start:
            j = d - next_i
            var_hi = max(j - next_i, var_hi)  # prune axis

    return dp[m - 1][n - 1], dp


# seq1 = "A" * 10
# seq2 = "A" * 10
# score_loop_band, dp_matrix_loop_band = needleman_wunsch_banded_loops(
#     seq1, seq2, bandwidth=2
# )

# for row in dp_matrix_loop_band:
#     print(row)


# Example usage
m, n = 7, 7  # dimensions
bandwidth = 1  # diagonal bandwidth
banded_matrix = create_banded_mask(m, n, bandwidth)
new_mat = [[0 for _ in range(n)] for _ in range(m)]
score, dp = x_drop(m - 1, n - 1, new_mat)
print("-------")
print_matrix(dp)
