import time


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Initialize the scoring matrix
    n, m = len(seq1), len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]

    start_time = time.time()

    # Initialize the first row and column
    for i in range(n + 1):
        score[i][0] = i * -gap
    for j in range(m + 1):
        score[0][j] = j * -gap

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (
                match if seq1[i - 1] == seq2[j - 1] else mismatch
            )
            delete = score[i - 1][j] + gap
            insert = score[i][j - 1] + gap
            score[i][j] = max(diag, delete, insert)

    end_time = time.time()
    print(f"Time taken to fill the scoring matrix: {end_time - start_time:.6f} seconds")

    # Traceback to get the alignment
    align1, align2 = "", ""
    i, j = n, m

    traceback_start_time = time.time()

    # Agg splitter max() or min() and returns each individual expression
    while i != 0 or j != 0:  # or between end points
        if i == 0:  # take base case and make if condition
            if score[i][j] == j * -gap:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1
        elif j == 0:  # take base case and make if condition
            if score[i][j] == i * -gap:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
        else:
            if score[i][j] == score[i - 1][j - 1] + (
                match if seq1[i - 1] == seq2[j - 1] else mismatch
            ):
                align1 = seq1[i - 1] + align1
                align2 = seq2[j - 1] + align2
                i -= 1
                j -= 1
            elif score[i][j] == score[i - 1][j] + gap:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
            elif score[i][j] == score[i][j - 1] + gap:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1

    traceback_end_time = time.time()
    print(
        f"Time taken for traceback: {traceback_end_time - traceback_start_time:.6f} seconds"
    )

    return align1, align2, score[n][m]


string_a = "A" * 3000
string_b = "B" * 3000
needleman_wunsch(string_a, string_b, match=1, mismatch=-1, gap=-2)
