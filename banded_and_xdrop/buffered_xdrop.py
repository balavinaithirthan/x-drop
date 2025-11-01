import matplotlib.pyplot as plt


def viz(mat, A, B):
    fig, ax = plt.subplots()
    ax.matshow(mat, cmap="RdYlGn")

    # Display the matrix values
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            c = mat[i][j]
            ax.text(j, i, str(c), va="center", ha="center")

    # Set tick labels to start at index 1
    ax.set_xticks(range(len(B) + 1))
    ax.set_yticks(range(len(A) + 1))

    ax.set_xticklabels([""] + list(B), fontsize=10)
    ax.set_yticklabels([""] + list(A), fontsize=10)

    plt.xlabel("B Sequence")
    plt.ylabel("A Sequence")
    plt.show()

def get_i_start(mat, d, i_start, i_end):
    i = i_start
    j = d - i
    while i < i_end and mat[i][j] == float("-inf"):
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
        explored[i][j] = 1
        if s < max_score - X_drop:
            mat[i][j] = float("-inf")
        # if (i == 4 and j == 1) or (i == 1 and j == 4):
        #     mat[i][j] = float('-inf')
        # if (i == 0 and j == 2) or (i == 2 and j == 0):
        #     mat[i][j] = 99


def modify_i_start(i_start, d_lo, ad):
    i_on_diag = (d_lo + ad) // 2
    i_start = max(i_start, i_on_diag)
    return i_start


def modify_i_end(i_end, d_hi, ad):
    i_on_diag = (d_hi + ad) // 2
    i_end = min(i_end, i_on_diag)
    return i_end

def nw(rows, cols):
    d_lo = 0
    d_hi = 0
    frontier_t_1 = [0 for _ in range(rows + cols)]
    frontier_t_2 = [0 for _ in range(rows + cols)]
    frontier_t_3 = [0 for _ in range(rows + cols)]
    for t in range(0, rows + cols + 1):
        for d in range(d_lo, d_hi + 1):
            frontier_t_1[d] = max(
                frontier_t_1[d], frontier_t_2[d - 1], frontier_t_3[d + 1]
            )
            markX(mat, t, i_start, i_end, x_thresh, max_score)
            i_start, d_lo = get_i_start(mat, ad, i_start, i_end)
            i_end, d_hi = get_i_end(mat, ad, i_start, i_end)
            lower_diag = max(lower_diag, d_lo)
            upper_diag = min(upper_diag, d_hi)
            i_start = modify_i_start(i_start, lower_diag, ad)
            i_end = modify_i_end(i_end, upper_diag, ad)
            d_lo = 
            d_hi = 

