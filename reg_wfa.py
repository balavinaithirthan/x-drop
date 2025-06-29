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
    max_e = m + n
    offset = max_e  # to map k=0 to center index

    # Allocate two wavefront buffers
    W = [[-1] * (2 * max_e + 1) for _ in range(2)]
    W[0][offset] = 0

    # Base case: extend matches from (0, 0)
    for k in range(0, 1):
        idx = offset + k
        v = W[0][idx] - k
        h = W[0][idx]
        while v < m and h < n and P[v] == T[h]:
            v += 1
            h += 1
            W[0][idx] += 1

    for e in range(1, max_e + 1):
        prev, curr = W[0], W[1]

        # Compute next wavefront at score e
        for k in range(-e, e + 1):
            idx = offset + k
            max_offset = -1
            if k + 1 <= e:
                max_offset = max(max_offset, prev[offset + (k + 1)])  # deletion
            if k - 1 >= -e:
                max_offset = max(max_offset, prev[offset + (k - 1)] + 1)  # insertion
            max_offset = max(max_offset, prev[idx] + 1)  # mismatch
            curr[idx] = max_offset

        # Extend matches
        for k in range(-e, e + 1):
            idx = offset + k
            v = curr[idx] - k
            h = curr[idx]
            while v < m and h < n and P[v] == T[h]:
                v += 1
                h += 1
                curr[idx] += 1

        # Check if goal diagonal has reached end of pattern
        if curr[offset + (m - n)] >= m:
            return e

        # Swap buffers
        W[0], W[1] = W[1], W[0]

    return -1  # fallback, should not occur


P = "GATTACA"
T = "GAATA"
print("Edit distance:", wfa_edit_distance(P, T))  # Output: 2
