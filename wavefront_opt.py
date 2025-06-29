import numpy as np
import matplotlib.pyplot as plt


# visualize the matrix using matddplotlib
def viz(mat, rowSeq, colSeq):
    fig, ax = plt.subplots()
    ax.matshow(mat, cmap="RdYlGn")

    # Display the matrix values
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            c = mat[i][j]
            ax.text(j, i, str(c), va="center", ha="center")

    # Set tick labels to start at index 1
    ax.set_xticks(range(len(colSeq) + 1))
    ax.set_yticks(range(len(rowSeq) + 1))

    ax.set_xticklabels([""] + list(colSeq), fontsize=10)
    ax.set_yticklabels([""] + list(rowSeq), fontsize=10)

    plt.xlabel("col Sequence")
    plt.ylabel("row Sequence")
    plt.show()


# print matrix function
def print_matrix(mat):
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            print(mat[i][j], end=" ")
        print()


def score(mat, A, B, i, j):
    # diag = mat[i - 1][j - 1] + int(A[i - 1] == B[j - 1])
    # left = mat[i][j - 1] - 2
    # up = mat[i - 1][j] - 2
    # return max(diag, left, up)
    return 10


def get_i_start(mat, d, i_start, i_end):
    i = i_start
    j = d - i
    while i < i_end and mat[i][j] == float("-inf"):
        ##explored[i][j] = 1
        i += 1
        j = d - i
    if mat[i_start][d - i_start] != float("-inf"):  # if no -inf yet then set to seen to
        d_lo = float("-inf")
    else:
        d_lo = i - j
        # -2 is optional, here we are saying make the
        # it looks a bit nicer because we are cutting off the ad wrong
    return i, d_lo


def get_i_end(mat, d, i_start, i_end):
    i = i_end - 1  # right now up to and not including!
    j = d - i
    while i > i_start and mat[i][j] == float("-inf"):
        ##explored[i][j] = 1
        i -= 1
        j = d - i
    if mat[i_end - 1][d - (i_end - 1)] != float(
        "-inf"
    ):  # if no -inf yet then set to seen to
        d_hi = float("inf")
    else:
        d_hi = i - j + 2
    return i + 1, d_hi


def markX(mat, d, i_start, i_end, X_drop, max_score):
    print("markX", d, i_start, i_end)
    for i in range(i_start, i_end):
        j = d - i
        s = mat[i][j]
        # explored[i][j] = 1
        if s < max_score - X_drop:
            mat[i][j] = float("-inf")  # we can just move this to the scoring loop below


def modify_i_start(i_start, d_lo, ad):
    # if d_lo < 0:
    #     return i_start
    i_on_diag = (d_lo + ad) // 2
    i_start = max(i_start, i_on_diag)
    return i_start


def modify_i_end(i_end, d_hi, ad):
    # if d_hi < 0:
    #     return i_end
    i_on_diag = (d_hi + ad) // 2
    i_end = min(i_end, i_on_diag)
    return i_end


def diagToOffsetIndex(diag: int, rows: int, cols: int) -> int:
    """
    Maps a diagonal index (where diag = row - col) to an offset array index
    in a grid of size rows x cols.

    The offset array has rows + cols - 1 entries.
    Diagonal 0 maps to index = number of rows - 1 if rows <= cols.
    """
    center = rows - 1 if rows <= cols else cols - 1
    return diag + center


# NO BASE CASE INCLUDED
# m x n matrix, m is rows, n is cols


"""
Some notes: 
- the diagonal index skips faster than each antidiagonal because of the way its structured
- one way to avoid this is to set an out of bounds for the offset array or out of bounds for (i,j)
- or if reach end, done
- that is why they have a while true loop and not for ad in range
"""


def nw4(rowSeq, colSeq, rows, cols, mat, x_thresh):
    m = rows
    n = cols
    rows = rows + 1
    cols = cols + 1
    i_start = 0
    i_end = 1
    diag_lo = 0
    diag_hi = 0
    offsets = [0 for _ in range(rows + cols)]
    while True:
        # i_start = max(max(1, ad - m), i_start)
        # i_end = min(min(ad, n) + 1, i_end)
        for diag in range(diag_lo, diag_hi + 1):
            offsetIdx = diagToOffsetIndex(diag, rows, cols)
            print("diag", diag, offsetIdx)
            i = offsets[offsetIdx]
            j = diag + i
            print("i", i, "j", j)
            s = score(mat, rowSeq, colSeq, i, j)
            mat[i][j] = s
            offsets[offsetIdx] = i + 1
            offsets[offsetIdx + 1] = i
            offsets[offsetIdx - 1] = i + 1
            if (i, j) == (rows - 1, cols - 1):
                print("done")
                viz(mat, rowSeq, colSeq)
                return mat
            # viz(mat, rowSeq, colSeq)
            print(offsets)
        diag_lo -= 1
        diag_hi += 1


def wavefront_edit_distance(seq1, seq2):
    n, m = len(seq1), len(seq2)
    max_score = n + m  # worst-case score (all inserts or deletes)

    # Each wavefront at score s maps diagonal k to offset (furthest i along diagonal)
    wavefronts = [{} for _ in range(max_score + 1)]

    # Base case: alignment starts at (0, 0) → diagonal k=0, offset i=0
    wavefronts[0][0] = 0

    for s in range(max_score + 1):
        wf = wavefronts[s]
        for k in wf:
            i = wf[k]
            j = i - k

            # extend matches as far as possible
            while i < n and j < m and seq1[i] == seq2[j]:
                i += 1
                j += 1
            wf[k] = i  # update furthest reaching point

            # check if reached bottom-right
            if i >= n and j >= m:
                return s

        # prepare next wavefront (s+1)
        if s + 1 <= max_score:
            next_wf = wavefronts[s + 1]
            for k in wf:
                i = wf[k]

                # Insertion → (i, j+1) → diagonal k-1
                next_wf[k - 1] = max(next_wf.get(k - 1, 0), i)

                # Deletion → (i+1, j) → diagonal k+1
                next_wf[k + 1] = max(next_wf.get(k + 1, 0), i + 1)

                # Mismatch → (i+1, j+1) → same diagonal
                next_wf[k] = max(next_wf.get(k, 0), i + 1)
    return -1  # should never reach here if sequences are finite


rowSeq = "ACGGGG"
colSeq = "ACGGGG"
rows = len(rowSeq)
cols = len(colSeq)
mat = np.zeros((rows + 1, cols + 1))
matp = mat[1:, 1:]
matp.fill(0.0)
# base case
for i in range(1, rows + 1):
    mat[i][0] = -i
for j in range(1, cols + 1):
    mat[0][j] = -j


# mat = nw4(rowSeq, colSeq, rows, cols, mat, 3)
wavefront_edit_distance(rowSeq, colSeq)
# print_matrix(mat)
