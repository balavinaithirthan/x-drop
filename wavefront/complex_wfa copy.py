from typing import List


# penalties
O = 10  # gap-open
E = 2  # gap-extend
X = 1  # mismatch
# match = 0
shift = 500


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


def wf_align(
    q: str,
    t: str,
) -> int:
    m, n = len(q), len(t)
    A_k = n - m
    A_offset = m

    M = I = D = [[0 for _ in range(1000)] for _ in range(1000)]
    M_hi = M_lo = I_hi = I_lo = D_hi = D_lo = [0 for _ in range(1000)]
    M_hi[0] = M_lo[0] = I_hi[0] = I_lo[0] = D_hi[0] = D_lo[0] = 0
    # initial wave-front (score 0)
    M[0][0] = 0
    s = 0

    # main loop
    while True:
        wf_extend(M[s], q, t, M_lo[s], M_hi[s])

        if M[s][A_k] > A_offset:  # reached (m,n)
            return s  # edit score found

        s += 1
        wf_next(M, I, D, M_lo, M_hi, I_lo, I_hi, D_lo, D_hi, s)


# -------------------------------------------------------------
#  Needleman–Wunsch / Gotoh alignment with affine gap penalties
#  M[v][h] – cost ending with a match / mismatch
#  I[v][h] – cost ending inside an insertion      (gap in seq1)
#  D[v][h] – cost ending inside a deletion        (gap in seq2)
# -------------------------------------------------------------
from math import inf


def align_affine(seq1: str, seq2: str) -> float:
    """
    Return (M, I, D, optimal_cost) for global alignment of seq1→seq2
    """
    o = O
    e = E
    mismatch_cost = X
    m, n = len(seq1), len(seq2)

    # --- allocate full DP tables ------------------------------
    M = [[inf] * (n + 1) for _ in range(m + 1)]
    I = [[inf] * (n + 1) for _ in range(m + 1)]
    D = [[inf] * (n + 1) for _ in range(m + 1)]
    M[0][0] = 0  # start

    # --- initialise first row (v = 0) ------------------------
    for h in range(1, n + 1):
        # opening a single insertion at h = 1, then extending
        I[0][h] = (o + e) if h == 1 else I[0][h - 1] + e
        M[0][h] = I[0][h]  # can only come from I

    # --- initialise first column (h = 0) --------------------
    for v in range(1, m + 1):
        D[v][0] = (o + e) if v == 1 else D[v - 1][0] + e
        M[v][0] = D[v][0]  # can only come from D

    # --- main DP -------------------------------------------
    for v in range(1, m + 1):
        for h in range(1, n + 1):
            I[v][h] = min(M[v][h - 1] + o + e, I[v][h - 1] + e)  # open  # extend

            D[v][h] = min(M[v - 1][h] + o + e, D[v - 1][h] + e)  # open  # extend

            sub = 0 if seq1[v - 1] == seq2[h - 1] else mismatch_cost
            M[v][h] = min(I[v][h], D[v][h], M[v - 1][h - 1] + sub)  # match / mismatch

    return M[m][n]  # optimal score in bottom-right


def run_known_tests():
    test_cases = [
        # (seq1, seq2, expected_score)
        ("", "", 0),
        ("A", "A", 0),
        ("A", "G", 1),  # mismatch
        ("A", "", 10),  # single gap open
        ("A", "AA", 12),  # open + extend
        ("ACGT", "ACGT", 0),  # perfect match
        ("ACGT", "AGGT", 1),  # 1 mismatch
        ("ACGT", "A-GT", 10),  # 1 deletion
        ("ACGT", "AGGTT", 13),  # mismatch + insertion
        ("GATTACA", "GCATGCU", 5),  # classic example
    ]

    passed = 0
    for idx, (s1, s2, expected) in enumerate(test_cases):
        affine_score = align_affine(s1, s2)
        # wfa_score = wf_align(s1, s2)
        wfa_score = 0

        if affine_score == expected == wfa_score:
            print(f"✓ Test {idx+1}: PASS")
            passed += 1
        else:
            print(f"✗ Test {idx+1}: FAIL")
            print(f"  seq1 = {s1}\n  seq2 = {s2}")
            print(
                f"  Expected = {expected}, Affine = {affine_score}, WFA = {wfa_score}"
            )

    print(f"\n{passed}/{len(test_cases)} known test cases passed.")


# Run the tests
if __name__ == "__main__":
    run_known_tests()
