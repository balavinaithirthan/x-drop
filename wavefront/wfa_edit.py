import difflib


def print_wavefronts(wavefronts):
    # labels columns by k
    # for i in range(-len(wavefronts[0]) // 2, len(wavefronts[0]) // 2 + 1):
    #     print(f"{i:2d}", end=" ")
    for i, wavefront in enumerate(wavefronts):
        print(f"Wavefront {i}: ", end="")
        for k in range(-len(wavefront) // 2, len(wavefront) // 2 + 1):
            print(f"{wavefront[k]:2d}", end=" ")
        print()
    print()


def print_single_wavefront(e, wavefront):
    print(f"Wavefront {e}: ", end="")
    # labels columns by k
    # for i in range(-len(wavefront) // 2, len(wavefront) // 2 + 1):
    #     print(f"{i:2d}", end=" ")
    print("Wavefront: ", end="")
    for k in range(len(wavefront)):
        print(f"{wavefront[k]:2d}", end=" ")
    print()
    print()


def wfa_edit_distance(P, T):
    m, n = len(P), len(T)
    max_e = m + n  # upper bound on edit distance

    # Initialize wavefront arrays: 2 rows (current and next), each with size 2*max_e + 1
    W = [
        [-1] * (2 * max_e + 1) for _ in range(max_e + 1)
    ]  # W[0] = current, W[1] = next
    offset = max_e  # shift so diagonal k=0 maps to index[offset]

    # Base case: offset 0 on diagonal 0 (k = 0)
    W[0][offset + 0] = 0
    extend(P, T, W[0], offset, 0)

    e = 0
    while True:
        if W[e][offset + (m - n)] >= m:
            return W, e

        compute_next(W[e], W[e + 1], e, offset)
        e += 1
        extend(P, T, W[e], offset, e)


def compute_next(W_e, W_e_1, e, offset):
    for k in range(-e, e + 2):
        idx = offset + k
        W_e_1[idx] = max(
            W_e[idx - 1] + 1, W_e[idx] + 1, W_e[idx + 1]
        )  # Deletion, Mismatch, Insertion


def extend(P, T, W_e, offset, e):
    for k in range(-e, e + 1):
        idx = offset + k
        v = W_e[idx] - k
        h = W_e[idx]
        while v < len(P) and h < len(T) and P[v] == T[h]:
            v += 1
            h += 1
            W_e[idx] += 1


def compute_W_2(P, T):
    m = len(P)
    n = len(T)
    numDiags = m + n + 1  # number of diagonals, corresponds to max_e
    max_e = n  # maximum possible edit distance
    diag_offset_wfa_arr = m  # offset is where the diagonal 0 is, ie the number of rows (negative diagonals)
    diag_end = n - m + diag_offset_wfa_arr
    W = [[-1] * (numDiags) for _ in range(max_e + 1)]
    W[0][diag_offset_wfa_arr] = 0  # W[0][0] = 0
    v = W[0][diag_offset_wfa_arr]  # make offset for 0,0 is diagonal 0
    h = W[0][diag_offset_wfa_arr]
    while v < len(P) and h < len(T) and P[v] == T[h]:
        v += 1
        h += 1
        W[0][diag_offset_wfa_arr + 0] += 1
    off = W[0][diag_end]
    if off == n:
        return W, 0
    for e in range(1, max_e + 1):
        # compute value
        for k in range(-e, e + 1):
            idx = diag_offset_wfa_arr + k
            W[e][idx] = max(
                W[e - 1][idx - 1] + 1, W[e - 1][idx] + 1, W[e - 1][idx + 1]
            )  # deletion, mismatch, insertion
        # extend W[e][k]
        for k in range(-e, e + 1):
            idx = diag_offset_wfa_arr + k
            v = W[e][idx] - k
            h = W[e][idx]

            while v < len(P) and h < len(T) and P[v] == T[h]:
                v += 1
                h += 1
                W[e][idx] += 1
            # check only diagonal 0
        off = W[e][diag_end]
        if off == n:
            return W, e

    return W, e  # fallback in case no match found


def gpt_edit_distance(P, T):
    """
    Compute edit distance between P and T using Wavefront Alignment
    with a 2D array W[e][k] where k is shifted by `offset`.
    """
    n, m = len(P), len(T)
    # maximum possible edit distance
    max_e = n + m
    # W[e][k+offset] = furthest h reached at edit-cost e on diagonal k
    W = [[-1] * (2 * max_e + 1) for _ in range(max_e + 1)]
    offset = max_e

    def extend_row(e):
        """Greedily advance matches on wavefront row W[e]."""
        for k in range(-e, e + 1):
            idx = k + offset
            h = W[e][idx]
            if h < 0:
                continue
            v = h - k
            # slide along diagonal while characters match
            while v < n and h < m and P[v] == T[h]:
                v += 1
                h += 1
            W[e][idx] = h

    # target diagonal where (v,h)=(n,m) satisfies h - v = (m-n)
    k_target = m - n

    # initialize e = 0
    W[0][offset] = 0
    extend_row(0)
    if W[0][k_target + offset] == m:
        return W, 0

    # grow wavefront
    for e in range(1, max_e + 1):
        row_prev = W[e - 1]
        row_cur = W[e]
        # computeNext()
        for k in range(-e, e + 1):
            idx = k + offset
            best = -1

            # Deletion: from k+1 at cost e-1
            h_del = row_prev[idx + 1]
            if h_del >= 0:
                best = max(best, h_del)

            # Mismatch/Substitution: from k at cost e-1
            h_sub = row_prev[idx]
            if h_sub >= 0:
                best = max(best, h_sub + 1)

            # Insertion: from k-1 at cost e-1
            h_ins = row_prev[idx - 1]
            if h_ins >= 0:
                best = max(best, h_ins + 1)

            row_cur[idx] = best

        # absorb all possible matches
        extend_row(e)

        # check for completion
        if W[e][k_target + offset] == m:
            return W, e

    # fallback (should never hit if max_e = n+m)
    return W, e


def true_edit_distance(P, T):
    m, n = len(P), len(T)

    # Create a (m+1) x (n+1) matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if P[i - 1] == T[j - 1]:
                cost = 0
            else:
                cost = 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,  # Deletion
                dp[i][j - 1] + 1,  # Insertion
                dp[i - 1][j - 1] + cost,  # Substitution
            )

    return dp[m][n]


def compute_W(P, T):
    """
    Compute the wavefront alignment matrix W and edit distance between strings P and T.
    Returns:
        W: 2D list representing the wavefront matrix.
        d: edit distance (int).
    """
    m, n = len(P), len(T)
    max_e = max(m, n)
    num_diags = 2 * max_e + 1  # number of diagonals
    offset = m
    diag_end = offset + (n - m)
    # initialize wavefront matrix with -1 (unreached)
    W = [[-1] * num_diags for _ in range(max_e + 1)]
    # initial match (wavefront e=0)
    v = h = 0
    W[0][offset] = 0
    while v < m and h < n and P[v] == T[h]:
        v += 1
        h += 1
        W[0][offset] += 1

    # check for perfect prefix match covering full alignment
    if W[0][diag_end] == n:
        return W, 0

    # iterative wavefront expansion
    for e in range(1, max_e + 1):
        k_min, k_max = -e, e
        # compute best predecessor
        for k in range(k_min, k_max + 1):
            idx = offset + k
            # best = -1
            # # mismatch (stay on same diagonal)
            # if 0 <= idx < num_diags and W[e - 1][idx] >= 0:
            #     best = W[e - 1][idx] + 1
            # # insertion (move from k-1)
            # if idx - 1 >= 0 and W[e - 1][idx - 1] >= 0 and W[e - 1][idx - 1] + 1 > best:
            #     best = W[e - 1][idx - 1] + 1
            # # deletion (move from k+1)
            # if (
            #     idx + 1 < num_diags
            #     and W[e - 1][idx + 1] >= 0
            #     and W[e - 1][idx + 1] > best
            # ):
            #     best = W[e - 1][idx + 1]

            # W[e][idx] = best
            W[e][idx] = max(
                W[e - 1][idx + 1],  # deletion
                W[e - 1][idx] + 1,  # insertion
                W[e - 1][idx - 1] + 1,  # mismatch
            )

            # extend matches from each wavefront cell
            # for k in range(k_min, k_max + 1):
            idx = offset + k
            val = W[e][idx]
            if val < 0:
                continue
            v = val - k
            h = val
            while v < m and h < n and P[v] == T[h]:
                v += 1
                h += 1
                W[e][idx] += 1

            # check if reached end of T at diagonal corresponding to (n-m)
            if 0 <= diag_end < num_diags and W[e][diag_end] == n:
                return W, e

    # fallback: maximum edits
    return W, e


if __name__ == "__main__":
    P = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    T = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
    true_e = true_edit_distance(P, T)
    print(f"True Edit Distance: {true_e}")
    W, e = compute_W(P, T)
    print("Trying Edit distance:", e)

    # compare W_e and W_e_2 values
    # for i in range(len(W)):
    #     for j in range(len(W[i])):
    #         if W[i][j] != W_2[i][j]:
    #             print(f"Mismatch at W[{i}][{j}]: {W[i][j]} vs {W_2[i][j]}")
    #         else:
    #             continue

    # print wavefronts

    # for i in range(len(W)):
    # print_single_wavefront(i, W[i])
