from typing import List, Dict


def print_matrix(M: List[List[int]]) -> None:
    for row in M:
        print(" ".join(f"{val:5}" for val in row))


def wfa_align(q: str, t: str, penalties: Dict[str, int]) -> int:
    x, o, e = penalties["x"], penalties["o"], penalties["e"]

    m, n = len(q), len(t)
    max_score: int = max(
        max(m, n) * x, max(m, n) * (o + e)
    )  # why is this the max_score
    num_diags = 2 * max_score + 1  # worst case
    offset = m

    INF = 10**9
    M = [[-INF] * (num_diags) for _ in range(max_score + 1)]
    I = [[-INF] * (num_diags) for _ in range(max_score + 1)]
    D = [[-INF] * (num_diags) for _ in range(max_score + 1)]

    M_lo = I_lo = D_lo = [0] * (max_score + 1)
    M_hi = I_hi = D_hi = [0] * (max_score + 1)
    A_k = offset + (m - n)
    A_offset = m

    ######
    M[0][offset] = 0
    v = h = 0
    while v >= 0 and h >= 0 and v < m and h < n and q[v] == t[h]:
        v += 1
        h += 1
        M[0][offset] += 1
    if M[0][A_k] == A_offset:
        return 0

    for score in range(1, max_score + 1):
        # hi ← max{𝑀̃ˢ⁻ˣʰⁱ, 𝑀̃ˢ⁻ᵒ⁻ᵉʰⁱ, Ĩˢ⁻ᵉʰⁱ, 𝒟̃ˢ⁻ᵉʰⁱ} + 1
        # lo ← min{𝑀̃ˢ⁻ˣˡᵒ, 𝑀̃ˢ⁻ᵒ⁻ᵉˡᵒ, Ĩˢ⁻ᵉˡᵒ, 𝒟̃ˢ⁻ᵉˡᵒ} − 1
        hi = (
            max(
                M_hi[score - x],
                M_hi[score - o - e],
                I_hi[score - e],
                D_hi[score - e],
            )
            + 1
        )
        lo = (
            min(
                M_lo[score - x],
                M_lo[score - o - e],
                I_lo[score - e],
                D_lo[score - e],
            )
            - 1
        )
        # lo =
        for k in range(lo, hi):
            idx = offset + k
            ins = max(M[score - o - e][idx - 1] + 1, I[score - e][idx - 1] + 1)
            delt = max(M[score - o - e][idx + 1], D[score - e][idx + 1])
            sub = M[score - x][idx] + 1

            I[score][idx] = ins
            D[score][idx] = delt
            M[score][idx] = max(sub, ins, delt, 0)

            val = M[score][idx]
            v = val - k
            h = val
            while v < m and h < n and q[v] == t[h]:
                v += 1
                h += 1
                M[score][idx] += 1

            if M[score][A_k] == A_offset:
                return score
        M_hi[score] = I_hi[score] = D_hi[score] = (
            max(M_hi[score - x], M_hi[score - o - e], I_hi[score - e], D_hi[score - e])
            + 1
        )
        M_lo[score] = I_lo[score] = D_lo[score] = (
            max(M_lo[score - x], M_lo[score - o - e], I_lo[score - e], D_lo[score - e])
            - 1
        )
    return 1000


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

    print_matrix(M)
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
        ("ACTGACTG", "", {"x": 1, "o": 2, "e": 1}),
    ]
    for q, t, penalties in test_cases:
        normal_score = normal_align(q, t, penalties)
        wfa_score = wfa_align(q, t, penalties)

        print(f"Testing sequences: {q} vs {t}")
        print(f"Penalties: {penalties}")
        print(f"Normal Align Score: {normal_score}")
        print(f"WFA Align Score: {wfa_score}")
        print("Match" if normal_score == wfa_score else "Mismatch")
        assert normal_score == wfa_score, f"Scores do not match for {q} vs {t}"
        print("-" * 40)


if __name__ == "__main__":
    test_alignments()
