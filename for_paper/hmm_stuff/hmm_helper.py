from math import log
from typing import Tuple, List

NEG_INF = float("-inf")

def safe_log(p: float) -> float:
    return NEG_INF if p <= 0.0 else log(p)

class PairHMMParams:
    """
    Minimal 3-state pair-HMM (Match M, Insertion I, Deletion D) with:
      - Uniform backgrounds for I/D emissions
      - Uniform per-letter equal/mismatch model for M emissions
      - Simple transition scheme with gap-open (M->I/D) and gap-extend (I->I, D->D)

    Start transitions come from a silent S to {M, I, D}.
    Stop is implicit (we just take max of M/I/D at (m,n)).
    """
    def __init__(
        self,
        alphabet: str = "ACGT",
        p_eq: float = 0.9,
        gap_open: float = 0.05,
        gap_extend: float = 0.4,
        start_open: float | None = None,
    ):
        assert 0.0 < p_eq < 1.0, "p_eq must be in (0,1)"
        if start_open is None:
            start_open = gap_open
        self.alphabet = alphabet
        self.A = len(alphabet)

        # Emission model
        # - M state: probability mass p_eq on equal pairs (spread uniformly),
        #            and (1 - p_eq) on mismatches (spread uniformly).
        self.eq_pair_log = safe_log(p_eq / self.A)
        self.mismatch_pair_log = safe_log((1.0 - p_eq) / (self.A * (self.A - 1)))
        # - I and D states emit a single symbol against a gap; use uniform background
        self.log_bg = safe_log(1.0 / self.A)

        # Transition model (logs)
        self.a = {}
        # Start
        self.a[("S", "M")] = safe_log(max(1e-12, 1.0 - 2.0 * start_open))
        self.a[("S", "I")] = safe_log(start_open)
        self.a[("S", "D")] = safe_log(start_open)
        # From M
        a_MI = gap_open
        a_MD = gap_open
        a_MM = max(1e-12, 1.0 - a_MI - a_MD)
        self.a[("M", "M")] = safe_log(a_MM)
        self.a[("M", "I")] = safe_log(a_MI)
        self.a[("M", "D")] = safe_log(a_MD)
        # From I
        a_II = gap_extend
        a_IM = max(1e-12, 1.0 - a_II)
        self.a[("I", "I")] = safe_log(a_II)
        self.a[("I", "M")] = safe_log(a_IM)
        # From D
        a_DD = gap_extend
        a_DM = max(1e-12, 1.0 - a_DD)
        self.a[("D", "D")] = safe_log(a_DD)
        self.a[("D", "M")] = safe_log(a_DM)

def log_emit_M(params: PairHMMParams, x: str, y: str) -> float:
    return params.eq_pair_log if x == y else params.mismatch_pair_log

def log_emit_I(params: PairHMMParams, y: str) -> float:
    return params.log_bg

def log_emit_D(params: PairHMMParams, x: str) -> float:
    return params.log_bg

# -------------------------
# Full O(m·n) Viterbi DP
# -------------------------
def viterbi_pair_hmm_full(seq1: str, seq2: str, params: PairHMMParams):
    """
    Standard Viterbi for a 3-state pair-HMM (M, I, D).
    Returns:
      end_score: best log-prob of aligning seq1..seq2
      path:      string over {M, X, I, D} (M=match, X=mismatch, I=ins in seq2, D=del from seq1)
      matrices:  (M, I, D) DP tables for inspection
    """
    m, n = len(seq1), len(seq2)
    if m == 0 and n == 0:
        return 0.0, "", ([[]], [[]], [[]])

    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    I = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    D = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    # First row (I only)
    if n >= 1:
        I[0][1] = params.a[("S", "I")] + log_emit_I(params, seq2[0])
        for j in range(2, n + 1):
            I[0][j] = (I[0][j - 1] + params.a[("I", "I")]) + log_emit_I(params, seq2[j - 1])

    # First column (D only)
    if m >= 1:
        D[1][0] = params.a[("S", "D")] + log_emit_D(params, seq1[0])
        for i in range(2, m + 1):
            D[i][0] = (D[i - 1][0] + params.a[("D", "D")]) + log_emit_D(params, seq1[i - 1])

    # First match (1,1) may start from S
    if m >= 1 and n >= 1:
        M[1][1] = params.a[("S", "M")] + log_emit_M(params, seq1[0], seq2[0])

    # Fill DP
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # M
            cand_M = [
                M[i - 1][j - 1] + params.a[("M", "M")],
                I[i - 1][j - 1] + params.a[("I", "M")],
                D[i - 1][j - 1] + params.a[("D", "M")],
            ]
            if i == 1 and j == 1:
                cand_M.append(0.0 + params.a[("S", "M")])
            M[i][j] = max(cand_M) + log_emit_M(params, seq1[i - 1], seq2[j - 1])

            # I
            I[i][j] = max(
                M[i][j - 1] + params.a[("M", "I")],
                I[i][j - 1] + params.a[("I", "I")],
            ) + log_emit_I(params, seq2[j - 1])

            # D
            D[i][j] = max(
                M[i - 1][j] + params.a[("M", "D")],
                D[i - 1][j] + params.a[("D", "D")],
            ) + log_emit_D(params, seq1[i - 1])

    end_score = max(M[m][n], I[m][n], D[m][n])

    # Traceback without storing full pointers (recompute local argmax)
    path = []
    i, j = m, n
    state = "M" if end_score == M[m][n] else ("I" if end_score == I[m][n] else "D")

    while i > 0 or j > 0:
        if state == "M":
            if i == 0 or j == 0:
                break
            s1 = M[i - 1][j - 1] + params.a[("M", "M")]
            s2 = I[i - 1][j - 1] + params.a[("I", "M")]
            s3 = D[i - 1][j - 1] + params.a[("D", "M")]
            s0 = 0.0 + params.a[("S", "M")] if (i == 1 and j == 1) else NEG_INF
            prev = max([(s1, "M"), (s2, "I"), (s3, "D"), (s0, "S")], key=lambda x: x[0])[1]
            path.append("M" if seq1[i - 1] == seq2[j - 1] else "X")
            i -= 1
            j -= 1
            if prev == "S":
                break
            state = prev

        elif state == "I":
            s1 = M[i][j - 1] + params.a[("M", "I")]
            s2 = I[i][j - 1] + params.a[("I", "I")]
            prev = max([(s1, "M"), (s2, "I")], key=lambda x: x[0])[1]
            path.append("I")
            j -= 1
            state = prev

        elif state == "D":
            s1 = M[i - 1][j] + params.a[("M", "D")]
            s2 = D[i - 1][j] + params.a[("D", "D")]
            prev = max([(s1, "M"), (s2, "D")], key=lambda x: x[0])[1]
            path.append("D")
            i -= 1
            state = prev

        else:
            raise RuntimeError("Invalid state in traceback")

    path.reverse()
    return end_score, "".join(path), (M, I, D)

# ----------------------------------------------------
# Wavefront-style Viterbi (adaptive band + z-drop)
# ----------------------------------------------------
def viterbi_pair_hmm_wavefront(
    seq1: str,
    seq2: str,
    params: PairHMMParams,
    band_half: int | None = None,
    band_max: int | None = None,
    z_drop: float = 50.0,
):
    """
    Adaptive-banded (wavefront-like) Viterbi in diagonal (k = i - j) space.

    Strategy:
      - Keep a moving diagonal band centered on the best-scoring diagonal seen so far
      - Only compute cells inside that band
      - Apply simple z-drop pruning: drop cells with score < (best - z_drop)

    With a very wide band and large z_drop, this matches the full DP.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 and n == 0:
        return 0.0, "", ([[]], [[]], [[]])

    if band_half is None:
        band_half = max(8, abs(m - n) + 2)
    if band_max is None:
        band_max = max(m, n)
    band_half = min(band_half, band_max)

    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    I = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    D = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    # Initialize first row/col (same as full DP so starting gaps are allowed)
    if n >= 1:
        I[0][1] = params.a[("S", "I")] + log_emit_I(params, seq2[0])
        for j in range(2, n + 1):
            I[0][j] = (I[0][j - 1] + params.a[("I", "I")]) + log_emit_I(params, seq2[j - 1])
    if m >= 1:
        D[1][0] = params.a[("S", "D")] + log_emit_D(params, seq1[0])
        for i in range(2, m + 1):
            D[i][0] = (D[i - 1][0] + params.a[("D", "D")]) + log_emit_D(params, seq1[i - 1])
    if m >= 1 and n >= 1:
        M[1][1] = params.a[("S", "M")] + log_emit_M(params, seq1[0], seq2[0])

    best_score = max(
        M[1][1] if m >= 1 and n >= 1 else NEG_INF,
        I[0][1] if n >= 1 else NEG_INF,
        D[1][0] if m >= 1 else NEG_INF,
        0.0,
    )
    best_diag = 0  # k = i - j

    for i in range(0, m + 1):
        # restrict j by current diagonal band
        j_min = max(0, i - (best_diag + band_half))
        j_max = min(n, i - (best_diag - band_half))
        if j_min > j_max:
            continue

        for j in range(j_min, j_max + 1):
            if i == 0 and j == 0:
                continue

            # M
            if i > 0 and j > 0:
                m_cand = [
                    M[i - 1][j - 1] + params.a[("M", "M")],
                    I[i - 1][j - 1] + params.a[("I", "M")],
                    D[i - 1][j - 1] + params.a[("D", "M")],
                ]
                if i == 1 and j == 1:
                    m_cand.append(0.0 + params.a[("S", "M")])
                M[i][j] = max(m_cand) + log_emit_M(params, seq1[i - 1], seq2[j - 1])

            # I
            if j > 0:
                I[i][j] = max(
                    M[i][j - 1] + params.a[("M", "I")],
                    I[i][j - 1] + params.a[("I", "I")],
                ) + log_emit_I(params, seq2[j - 1])

            # D
            if i > 0:
                D[i][j] = max(
                    M[i - 1][j] + params.a[("M", "D")],
                    D[i - 1][j] + params.a[("D", "D")],
                ) + log_emit_D(params, seq1[i - 1])

            # Update best & z-drop prune
            cell_best = max(M[i][j], I[i][j], D[i][j])
            if cell_best > best_score:
                best_score = cell_best
                best_diag = i - j
            if best_score - cell_best > z_drop:
                M[i][j] = I[i][j] = D[i][j] = NEG_INF

    end_score = max(M[m][n], I[m][n], D[m][n])

    # Traceback (same logic as full DP)
    path = []
    i, j = m, n
    state = "M" if end_score == M[m][n] else ("I" if end_score == I[m][n] else "D")

    while i > 0 or j > 0:
        if state == "M":
            if i == 0 or j == 0:
                break
            s1 = M[i - 1][j - 1] + params.a[("M", "M")]
            s2 = I[i - 1][j - 1] + params.a[("I", "M")]
            s3 = D[i - 1][j - 1] + params.a[("D", "M")]
            s0 = 0.0 + params.a[("S", "M")] if (i == 1 and j == 1) else NEG_INF
            prev = max([(s1, "M"), (s2, "I"), (s3, "D"), (s0, "S")], key=lambda x: x[0])[1]
            path.append("M" if seq1[i - 1] == seq2[j - 1] else "X")
            i -= 1
            j -= 1
            if prev == "S":
                break
            state = prev

        elif state == "I":
            s1 = M[i][j - 1] + params.a[("M", "I")]
            s2 = I[i][j - 1] + params.a[("I", "I")]
            prev = max([(s1, "M"), (s2, "I")], key=lambda x: x[0])[1]
            path.append("I")
            j -= 1
            state = prev

        elif state == "D":
            s1 = M[i - 1][j] + params.a[("M", "D")]
            s2 = D[i - 1][j] + params.a[("D", "D")]
            prev = max([(s1, "M"), (s2, "D")], key=lambda x: x[0])[1]
            path.append("D")
            i -= 1
            state = prev

        else:
            raise RuntimeError("Invalid state in traceback")

    path.reverse()
    return end_score, "".join(path), (M, I, D)

# -------------------------
# Pretty-print alignment
# -------------------------
def render_alignment(seq1: str, seq2: str, path: str) -> Tuple[str, str, str]:
    """Return (aligned_seq1, midline, aligned_seq2)."""
    i = j = 0
    a1, a2, mid = [], [], []
    for op in path:
        if op in ("M", "X"):
            a1.append(seq1[i]); a2.append(seq2[j])
            mid.append("|" if seq1[i] == seq2[j] else ".")
            i += 1; j += 1
        elif op == "I":
            a1.append("-"); a2.append(seq2[j]); mid.append(" ")
            j += 1
        elif op == "D":
            a1.append(seq1[i]); a2.append("-"); mid.append(" ")
            i += 1
        else:
            raise ValueError(f"Unknown op {op}")
    return "".join(a1), "".join(mid), "".join(a2)

# -------------------------
# Quick tests / examples
# -------------------------
if __name__ == "__main__":
    params = PairHMMParams(
        alphabet="ACGT",
        p_eq=0.9,
        gap_open=0.05,
        gap_extend=0.4,
    )

    tests = [
        ("", ""),
        ("A", ""),
        ("", "A"),
        ("A", "A"),
        ("ACGT", "ACGT"),
        ("ACGT", "AGGT"),
        ("GATTACA", "GCATGCU"),
        ("ACGTACGT", ""),
        ("", "ACGTACGT"),
        ("AAAA", "TTTT"),
    ]

    for s1, s2 in tests:
        s_full, path_full, _ = viterbi_pair_hmm_full(s1, s2, params)
        s_wf, path_wf, _ = viterbi_pair_hmm_wavefront(s1, s2, params, band_half=50, z_drop=1e9)

        ok = abs(s_full - s_wf) < 1e-9
        print(f"{s1!r} vs {s2!r}:")
        print("  full      :", s_full, path_full)
        print("  wavefront :", s_wf, path_wf, "  (match)" if ok else "  (!! differs)")
        if s1 and s2:
            a1, mid, a2 = render_alignment(s1, s2, path_full)
            print("  aln (full):")
            print("   ", a1)
            print("   ", mid)
            print("   ", a2)
        print()
