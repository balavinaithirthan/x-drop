from typing import List, Dict

# NEG_INF = -(10**9)  # “does not exist”
# MAX_SCORE = 1000  # crude upper bound; raise if needed


# # ----------------------------------------------------------------------
# # WF_NEXT  –  build wave-front for score s
# # ----------------------------------------------------------------------
# def wf_next(
#     M: List[List[int]],
#     I: List[List[int]],
#     D: List[List[int]],
#     M_lo: List[int],
#     M_hi: List[int],
#     I_lo: List[int],
#     I_hi: List[int],
#     D_lo: List[int],
#     D_hi: List[int],
#     s: int,
#     X: int,
#     O: int,
#     E: int,
#     shift: int,
# ) -> None:
#     # 1. band limits for new score
#     hi = max(M_hi[s - X], M_hi[s - O - E], I_hi[s - E], D_hi[s - E]) + 1
#     lo = min(M_lo[s - X], M_lo[s - O - E], I_lo[s - E], D_lo[s - E]) - 1

#     # 3. populate every diagonal in the new band
#     for k in range(lo, hi + 1):  # inclusive!
#         ins = max(M[s - O - E][k - 1] + 1, I[s - E][k - 1] + 1)
#         delt = max(M[s - O - E][k + 1], D[s - E][k + 1])
#         sub = M[s - X][k] + 1

#         I[s][k + shift] = ins
#         D[s][k + shift] = delt
#         M[s][k + shift] = max(sub, ins, delt)  # Eq. 3

#     # M̃_hi_s = Ĩ_hi_s = D̃_hi_s = max{M̃_hi_{s−x}, M̃_hi_{s−o−e}, Ĩ_hi_{s−e}, D̃_hi_{s−e}} + 1
#     # M̃_lo_s = Ĩ_lo_s = D̃_lo_s = max{M̃_lo_{s−x}, M̃_lo_{s−o−e}, Ĩ_{s−e}, D̃_{s−e}} − 1
#     M_hi[s] = I_hi[s] = D_hi[s] = (
#         max(M_hi[s - X], M_hi[s - O - E], I_hi[s - E], D_hi[s - E]) + 1
#     )
#     M_lo[s] = I_lo[s] = D_lo[s] = (
#         max(M_lo[s - X], M_lo[s - O - E], I_lo[s - E], D_lo[s - E]) - 1
#     )


# def wf_extend(
#     M: List[int],
#     q: str,
#     t: str,
#     M_lo: int,
#     M_hi: int,
# ) -> None:
#     m, n = len(q), len(t)
#     for k in range(M_lo, M_hi):
#         v = M[k] - k  # row index
#         h = M[k]  # col index
#         while v < m and h < n and q[v] == t[h]:
#             M[k] += 1
#             v += 1
#             h += 1


# def wf_align(
#     q: str,
#     t: str,
#     x: int,
#     o: int,
#     e: int,
# ) -> int:
#     shift = max(len(q), len(t))  # to avoid negative indices
#     m, n = len(q), len(t)
#     A_k = n - m
#     A_offset = m

#     M = I = D = [[0 for _ in range(1000)] for _ in range(1000)]
#     M_hi = M_lo = I_hi = I_lo = D_hi = D_lo = [0 for _ in range(1000)]
#     M_hi[0] = M_lo[0] = I_hi[0] = I_lo[0] = D_hi[0] = D_lo[0] = 0
#     # initial wave-front (score 0)
#     M[0][0] = 0
#     s = 0

#     # main loop
#     while True:
#         wf_extend(M[s], q, t, M_lo[s], M_hi[s])

#         if M[s][A_k] > A_offset:  # reached (m,n)
#             return s  # edit score found

#         s += 1
#         wf_next(M, I, D, M_lo, M_hi, I_lo, I_hi, D_lo, D_hi, s, x, o, e, shift)


def wfa_align(q: str, t: str, penalties: Dict[str, int]) -> int:
    if q == t:
        djiwq = 0

    x, o, e = penalties["x"], penalties["o"], penalties["e"]
    m, n = len(q), len(t)
    INF = 10**9
    M = [[-INF] * (1000) for _ in range(1000)]
    I = [[-INF] * (1000) for _ in range(1000)]
    D = [[-INF] * (1000) for _ in range(1000)]
    offset = 150
    A_k = m - n + offset  # where the diagonal 0 is
    A_offset = m

    min_score = 0
    max_score = 1000
    min_diag = -50
    max_diag = 50

    ######
    M[0][offset] = 0
    v = h = 0
    while v >= 0 and h >= 0 and v < m and h < n and q[v] == t[h]:
        v += 1
        h += 1
        M[0][offset] += 1
    if M[0][A_k] == A_offset:
        return 0

    #######

    for score in range(1, 1000):
        for k in range(-50, 50):
            idx = offset + k
            ins = max(M[score - o - e][idx - 1] + 1, I[score - e][idx - 1] + 1)
            delt = max(M[score - o - e][idx + 1], D[score - e][idx + 1])
            sub = M[score - x][idx] + 1

            I[score][idx] = ins
            D[score][idx] = delt
            M[score][idx] = max(sub, ins, delt)

            val = M[score][idx]
            v = val - k
            h = val
            while v >= 0 and h >= 0 and v < m and h < n and q[v] == t[h]:
                v += 1
                h += 1
                M[score][idx] += 1

            if M[score][A_k] == A_offset:
                return score
    return 1000  # fallback: maximum edits


# ---------------------------------------------------------------------
#  Classic DP with affine gaps (min‑score Needleman–Wunsch + gaps)
# ---------------------------------------------------------------------
def normal_align(q: str, t: str, penalties: Dict[str, int]) -> int:
    """
    Traditional dynamic‑programming alignment with gap‑open (o) and
    gap‑extend (e) penalties and mismatch penalty x (match = 0).
    Returns M[n][m] – the final score in the M matrix.
    """
    x, o, e = penalties["x"], penalties["o"], penalties["e"]
    n, m = len(q), len(t)
    INF = 10**9

    # allocate (n+1) x (m+1) matrices
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    I = [[INF] * (m + 1) for _ in range(n + 1)]
    D = [[INF] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0

    # initialise first row (v = 0)
    for h in range(1, m + 1):
        I[0][h] = (o + e) + (h - 1) * e
        M[0][h] = I[0][h]

    # initialise first column (h = 0)
    for v in range(1, n + 1):
        D[v][0] = (o + e) + (v - 1) * e
        M[v][0] = D[v][0]

    # fill rest
    for v in range(1, n + 1):
        for h in range(1, m + 1):
            I[v][h] = min(M[v][h - 1] + o + e, I[v][h - 1] + e)
            D[v][h] = min(M[v - 1][h] + o + e, D[v - 1][h] + e)
            sub = M[v - 1][h - 1] + (0 if q[v - 1] == t[h - 1] else x)
            M[v][h] = min(I[v][h], D[v][h], sub)

    return M[n][m]


def test_alignments():
    test_cases = [
        ("ACGT", "ACGT", {"x": 1, "o": 2, "e": 1}),
        ("ACGT", "TGCA", {"x": 2, "o": 3, "e": 2}),
        ("AAAA", "TTTT", {"x": 1, "o": 1, "e": 1}),
        ("ACGTACGT", "ACGTACGT", {"x": 0, "o": 5, "e": 2}),
        ("GATTACA", "GATTACA", {"x": 1, "o": 2, "e": 1}),
        ("GATTACA", "GATTTCA", {"x": 2, "o": 3, "e": 2}),
        ("ACTG", "ACCG", {"x": 1, "o": 2, "e": 1}),
        ("ACTG", "ACT", {"x": 1, "o": 2, "e": 1}),
        ("ACTGACTG", "ACTG", {"x": 1, "o": 2, "e": 1}),
        ("", "", {"x": 1, "o": 2, "e": 1}),
        ("A", "T", {"x": 1, "o": 2, "e": 1}),
        ("ACGT", "ACGTT", {"x": 1, "o": 2, "e": 1}),
        ("ACGT", "AC", {"x": 1, "o": 2, "e": 1}),
    ]
    for q, t, penalties in test_cases:
        normal_score = normal_align(q, t, penalties)
        wfa_score = wfa_align(q, t, penalties)

        print(f"Testing sequences: {q} vs {t}")
        print(f"Penalties: {penalties}")
        print(f"Normal Align Score: {normal_score}")
        print(f"WFA Align Score: {wfa_score}")
        print("Match" if normal_score == wfa_score else "Mismatch")
        print("-" * 40)


if __name__ == "__main__":
    test_alignments()
