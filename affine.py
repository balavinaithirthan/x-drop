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


def score(mat, A, B, i, j):
    diag = mat[i - 1][j - 1] + int(A[i - 1] == B[j - 1])
    left = mat[i][j - 1] - 2
    up = mat[i - 1][j] - 2
    return max(diag, left, up)


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


# if (i == 4 and j == 1) or (i == 1 and j == 4):
#     mat[i][j] = float('-inf')
# if (i == 0 and j == 2) or (i == 2 and j == 0):
#     mat[i][j] = 99


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


# NO BASE CASE INCLUDED
# m x n matrix, m is rows, n is cols
def nw4(rowSeq, colSeq, rows, cols, mat, x_thresh):
    rows = rows + 1
    cols = cols + 1
    i_start = 1
    i_end = 1
    max_score = float("-inf")
    # TODO: figure out how to incorporate the base case into the compute and maybe start from 0?
    lower_diag = float("-inf")
    upper_diag = float("inf")
    # this has to be the long side
    for ad in range(2, rows + cols - 1):
        # change to 0 and i_start = 0, i_end = 0 if base case incorporated
        max_i_end = rows
        i_end = min(
            i_end + 1, max_i_end
        )  # keep increasing i by 1 until it reaches the the lowest point (rows)
        max_i_start = ad - cols + 1
        i_start = max(
            i_start, max_i_start
        )  # stay at 0 until start traversing past half way point of ad matrix
        # i_start = min(i_start, rows)  # rows is the max, don't need bc of stopping condition
        if i_start >= i_end:
            break
        for i in range(i_start, i_end):
            j = ad - i
            s = score(mat, rowSeq, colSeq, i, j)
            mat[i][j] = s
            # explored[i][j] = 1
            max_score = max(mat[i][j], max_score)
        markX(mat, ad, i_start, i_end, x_thresh, max_score)
        i_start, d_lo = get_i_start(mat, ad, i_start, i_end)
        i_end, d_hi = get_i_end(mat, ad, i_start, i_end)
        lower_diag = max(lower_diag, d_lo)
        upper_diag = min(upper_diag, d_hi)
        i_start = modify_i_start(i_start, lower_diag, ad)
        i_end = modify_i_end(i_end, upper_diag, ad)
        print(i_start, i_end)
    viz(mat, rowSeq, colSeq)
    return mat


def init_matrix(rows, cols):
    mat = np.zeros((rows + 1, cols + 1))
    matp = mat[1:, 1:]
    matp.fill(0.0)
    # base case
    for i in range(1, rows + 1):
        mat[i][0] = -i
    for j in range(1, cols + 1):
        mat[0][j] = -j
    return mat


rowSeq = "ACGGGGGGG"
colSeq = "ACGGGG"
rows = len(rowSeq)
cols = len(colSeq)
M = init_matrix(len(rowSeq), len(colSeq))
Ix = init_matrix(len(rowSeq), len(colSeq))
Iy = init_matrix(len(rowSeq), len(colSeq))


# viz(mat, rowSeq, colSeq)

mat = nw4(rowSeq, colSeq, rows, cols, M, 1)
# print_matrix(mat)


# three matrices
# MATRIX 1: Ix, M(i-1, j-1), Iy
# MATRIX 2: Iy
# MATRIX 3: M
