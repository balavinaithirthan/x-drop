from wfa_edit import compute_W, wfa_edit_distance, gpt_edit_distance, true_edit_distance
from Wavefront_Aligner_Class import WavefrontAligner


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

    # try:
    #     W2, dist2 = wfa_edit_distance(P, T)
    #     print(
    #         f"{'✅' if dist2 == true_dist else '❌'} wfa_edit_distance distance: {dist2}"
    #     )
    # except Exception as e:
    #     print(f"❌ wfa_edit_distance crashed: {e}")

    # try:
    #     W3, dist3 = gpt_edit_distance(P, T)
    #     print(
    #         f"{'✅' if dist3 == true_dist else '❌'} gpt_edit_distance distance: {dist3}"
    #     )
    # except Exception as e:
    #     print(f"❌ gpt_edit_distance crashed: {e}")


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
    ]
    test_cases += [
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


if __name__ == "__main__":
    run_all_validations()
