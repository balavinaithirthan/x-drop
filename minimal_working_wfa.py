def compute_W(P, T):
    m, n = len(P), len(T)
    max_e = max(m, n)  # the maximum possible edit distance/score
    num_diags = m + n + 1  # number of diagonals
    # this is the minimum we need, but if we just make this num_diags = 2 * max_e + 1, don't have to bounds check (look at case "ACGTACGT" vs "")
    num_diags = 2 * max_e + 1
    offset = m  # offset off main diagonal, ie diagonal 0 should be at idx offset in the wavefront
    diag_end = offset + (n - m)  # location of the diagonal the hits the (m,n) point
    W = [[-1] * num_diags for _ in range(max_e + 1)]  # score by number of diagonals
    v = h = 0
    W[0][offset] = 0
    while v < m and h < n and P[v] == T[h]:
        v += 1
        h += 1
        W[0][offset] += 1

    # check for perfect prefix match covering full alignment
    if W[0][diag_end] == n:
        return W, 0
    # wavefront expansion
    for e in range(1, max_e + 1):
        k_min, k_max = -e, e
        # compute best predecessor
        for k in range(k_min, k_max + 1):
            idx = offset + k
            W[e][idx] = max(
                W[e - 1][idx + 1],  # deletion
                W[e - 1][idx] + 1,  # insertion
                W[e - 1][idx - 1] + 1,  # mismatch
            )

            # extend matches from each wavefront cell
            # extend encodes the +0 for the matches
            val = W[e][idx]
            j = val - k  # j - i + i, so should be be k + val?
            i = val  # i
            while i < m and j < n and P[i] == T[j]:
                i += 1
                j += 1
                W[e][idx] += 1

            # check if reached end of T at diagonal corresponding to (n-m)
            if W[e][diag_end] == n:
                return W, e
    # fallback: maximum edits
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


def validate_edit_distances(P, T, test_label=""):
    """Run and compare edit distance results for a single pair of strings."""
    print(f"\n🧪 {test_label}")
    print(f"Pattern (P): {P}")
    print(f"Text    (T): {T}")

    true_dist = true_edit_distance(P, T)
    print(f"✅ True Edit Distance: {true_dist}")

    try:
        W1, dist1 = compute_W(P, T)
        print(f"{'✅' if dist1 == true_dist else '❌'} compute_W distance: {dist1}")
    except Exception as e:
        print(f"❌ compute_W crashed: {e}")


def run_all_validations():
    test_cases = [
        ("ACGT", "ACGT", "Exact match (tiny)"),
        ("ACGT", "AGCT", "1 substitution"),
        ("ACGT", "ACGTTT", "Insertion at end"),
        ("ACGT", "TTACGT", "Insertion at start"),
        ("", "ACGTACGT", "Empty vs full"),
        ("ACGTACGT", "", "Full vs empty"),
        ("ACGTACGT", "ACGT", "Truncation"),
        ("ACGTACGT", "AGGTTCGT", "Multiple substitutions"),
        ("ACGT" * 50, "ACGT" * 50, "Large identical strings"),
        ("ACGT" * 50, "TGCA" * 50, "Large reversed bases"),
        ("A" * 100, "A" * 99 + "T", "Single end substitution in long string"),
        ("A", "A", "Single character match"),  # Expected: 0
        ("A", "G", "Single character substitution"),  # Expected: 1
        ("", "", "Empty strings"),  # Expected: 0
        ("ACTG", "TGCA", "All substitutions"),  # Expected: 4
        ("ACGT", "ACGTACGT", "Exact prefix match"),  # Expected: 4
        ("ACGTACGT", "ACGT", "Exact suffix match removed"),  # Expected: 4
        ("GATTACA", "GCATGCU", "Classic Levenshtein example"),  # Expected: 4
        ("AAAGGG", "AAAGG", "One deletion at end"),  # Expected: 1
        ("AAAGGG", "AAGGG", "One deletion at start"),  # Expected: 1
        ("AAAGGG", "AAAAGGG", "One insertion at start"),  # Expected: 1
        ("AAAGGG", "AAAGGGG", "One insertion at end"),  # Expected: 1
        ("AGCTTAGCTA", "AGCTAGCTA", "Single deletion in middle"),  # Expected: 1
        ("ACGT" * 25, "ACGT" * 25 + "A", "Long with 1 extra char"),  # Expected: 1
        ("A" * 500, "A" * 499 + "C", "Long string single mismatch"),  # Expected: 1
        ("ACGT" * 30, "ACG" * 40, "Pattern length mismatch"),  # Expected: high
        ("ACGTACGT", "TGCATGCA", "All positions shifted"),  # Expected: 8
        ("XYZ", "xyz", "Case sensitivity"),  # Expected: 3
        ("ACGT", "ACXT", "One internal mismatch"),  # Expected: 1
        ("AGCTAGCTAG", "AGCTTGCTAG", "Substitution in the middle"),  # Expected: 1
    ]

    for P, T, label in test_cases:
        validate_edit_distances(P, T, test_label=label)


# if __name__ == "__main__":
#     run_all_validations()


if __name__ == "__main__":
    P = ""
    T = "ACGTACGT"
    true_e = true_edit_distance(P, T)
    print(f"True Edit Distance: {true_e}")
    W, e = compute_W(P, T)
    print("Trying Edit distance:", e)
