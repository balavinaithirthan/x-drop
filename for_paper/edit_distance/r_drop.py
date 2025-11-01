import math

def r_drop(query, target, Rthresh=4, match=1, mismatch=-1, gap=-1):
    """
    Repeat-aware drop (R-drop): prunes repetitive axes dynamically.
    Simple run-length heuristic version.
    """

    m, n = len(query), len(target)
    dp = [[-math.inf] * (n + 1) for _ in range(m + 1)]
    runlen = [[0] * (n + 1) for _ in range(m + 1)]

    dp[0][0] = 0
    var_lo, var_hi = math.inf, -math.inf
    best_score = 0

    for d in range(0, m + n + 1):  # anti-diagonals
        i_start = max(1, d - n)
        i_end   = min(m, d)

        for i in range(i_start, i_end + 1):
            j = d - i
            if j < 1 or j > n:
                continue

            # Compute local run-length score (repeat proxy)
            if query[i-1] == query[i-2] if i > 1 else False:
                runlen[i][j] = runlen[i-1][j]
            elif target[j-1] == target[j-2] if j > 1 else False:
                runlen[i][j] = runlen[i][j-1]
            else:
                runlen[i][j] = 0

            # Apply recurrence with repeat penalty
            score = max(
                dp[i-1][j-1] + (match if query[i-1] == target[j-1] else mismatch),
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

            # If too repetitive, prune
            if runlen[i][j] > Rthresh:
                dp[i][j] = -math.inf
                continue

            dp[i][j] = score
            best_score = max(best_score, score)

        # Update axis limits based on active cells
        active = [(d - i, i) for i in range(i_start, i_end + 1) if dp[i][d - i] != -math.inf]
        if active:
            var_lo = min(var_lo, active[0][0] - active[0][1])
            var_hi = max(var_hi, active[-1][0] - active[-1][1])

    return dp, dp[m][n], var_lo, var_hi, best_score
