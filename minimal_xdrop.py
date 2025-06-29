import numpy as np
import matplotlib.pyplot as plt


# visualize the matrix using matplotlib
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


def score(mat, A, B, i, j, gap=1):
    if j == 0:
        return i * -gap
    elif i == 0:
        return j * -gap
    else:
        diag = mat[i - 1][j - 1] + int(A[i - 1] == B[j - 1])
        left = mat[i][j - 1] - 2
        up = mat[i - 1][j] - 2
        return max(diag, left, up)


def markX(mat, d, i_start, i_end, X_drop, max_score):
    # print("markX", d, i_start, i_end)
    for i in range(i_start, i_end):
        j = d - i
        s = mat[i][j]
        # explored[i][j] = 1
        if s < max_score - X_drop:
            mat[i][j] = float("-inf")  # we can just move this to the scoring loop below


def print_antidiagonal(matrix, ad):
    rows = len(matrix)
    cols = len(matrix[0])
    for i in range(rows):
        j = ad - i
        if 0 <= j < cols:
            print(matrix[i][j])


from math import inf


# NO BASE CASE INCLUDED
# m x n matrix, m is rows, n is cols
def nw4(rowSeq, colSeq, rows, cols, mat, x_thresh):
    rows = rows + 1  # 7 x 7, for len 6 str
    cols = cols + 1  # rows from 0...6 and cols from 0...6
    max_score = float("-inf")
    var_lo = -inf
    var_hi = inf
    # prune axis is different from the actual i & j coordinates!
    for ad in range(0, rows + cols + 1):  # there are rows + cols + 1 antidiagonals
        i_start = max(0, ad - cols + 1)  # activate after middle ad
        i_start = max(
            i_start, (ad - var_hi) // 2
        )  # lower is better prune but that means higher i
        i_end = min(rows - 1, ad)  # activate until middle ad, then do row -1
        i_end = min(
            i_end, (ad - var_lo) // 2
        )  # higher is better prune but that means lower i_end
        if i_start >= i_end + 1:
            return mat
        print(i_start, i_end)
        for i in range(i_start, i_end + 1):
            j = ad - i
            s = score(mat, rowSeq, colSeq, i, j)
            mat[i][j] = s
            max_score = max(mat[i][j], max_score)
        markX(
            mat, ad, i_start, i_end + 1, x_thresh, max_score
        )  # marks -inf along the antidiagonal based on the x-drop threshold
        next_i = 0
        for i in range(i_end, i_start - 1, -1):
            j = ad - i
            if mat[i][j] != -inf:
                next_i = i
                break
        if next_i != i_end:
            j = ad - next_i
            var_lo = max(j - next_i, var_lo)  # prune axis, higher is better prune
            print("prune lo at", var_lo)

        next_i = 0
        for i in range(i_start, i_end + 1):
            j = ad - i
            if mat[i][j] != -inf:
                next_i = i
                break
        if next_i != i_start:
            j = ad - next_i
            var_hi = min(j - next_i, var_hi)  # prune axis, lower is better prune
            print("prune hi at", var_hi)
    return mat


rowSeq = "GGGGGGAAAAAA"
colSeq = "AAAAATAAAAAA"
rows = len(rowSeq)
cols = len(colSeq)
mat = np.zeros((rows + 1, cols + 1))
matp = mat[1:, 1:]
matp.fill(0.0)
# base case
# for i in range(1, rows + 1):
#     mat[i][0] = -i
# for j in range(1, cols + 1):
#     mat[0][j] = -j


mat = nw4(rowSeq, colSeq, rows, cols, mat, 3)
viz(mat, rowSeq, colSeq)

# print_matrix(mat)
