import numpy as np
import matplotlib.pyplot as plt

A = "ACGGGG"
B = "ATTCGG"
ping = 1


def initialize(A, B, fill_val):
    N = len(A)
    mat = np.zeros((N+ping, N+ping))
    matp = mat[1:,1:]
    matp.fill(fill_val)
    #base case
    for i in range(1, len(A)+ping):
        mat[i][0] = -i
    for j in range(1, len(B)+ping):
        mat[0][j] = -j  
    return mat

explored = initialize(A,B, 2)


#visualize the matrix using matplotlib
def viz(mat, A, B):
    fig, ax = plt.subplots()
    ax.matshow(mat, cmap='RdYlGn')
    
    # Display the matrix values
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            c = mat[i][j]
            ax.text(j, i, str(c), va='center', ha='center')

    # Set tick labels to start at index 1
    ax.set_xticks(range(len(B) + 1))
    ax.set_yticks(range(len(A) + 1))
    
    ax.set_xticklabels([''] + list(B), fontsize=10)
    ax.set_yticklabels([''] + list(A), fontsize=10)
    
    plt.xlabel('B Sequence')
    plt.ylabel('A Sequence')
    plt.show()


#print matrix function
def print_matrix(mat):
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            print(mat[i][j], end = " ")
        print()

def score(mat, A, B, i,j):
    diag = mat[i-1][j-1] + int(A[i-1] == B[j-1])
    left = mat[i][j-1] - 2
    up = mat[i-1][j] - 2
    return  max(diag, left, up)


def get_i_start(mat,d, i_start, i_end):
    i = i_start
    j = d - i
    while i < i_end and mat[i][j] == float('-inf'):
        explored[i][j] = 1
        i += 1
        j = d - i
    if mat[i_start][d-i_start] != float('-inf'): # if no -inf yet then set to seen to 
        d_lo= float('-inf')
    else:
        d_lo = i - j
        # -2 is optional, here we are saying make the
        # it looks a bit nicer because we are cutting off the ad wrong
    return i, d_lo

def get_i_end(mat,d, i_start, i_end):
    i = i_end-1 # right now up to and not including!
    j = d - i
    while i > i_start and mat[i][j] == float('-inf'):
        explored[i][j] = 1
        i -= 1
        j = d - i
    if mat[i_end - 1][d - (i_end - 1)] != float('-inf'): # if no -inf yet then set to seen to 
        d_hi = float('inf')
    else:
        d_hi = i - j + 2
    return i + 1, d_hi

def markX(mat,d, i_start, i_end, X_drop, max_score):
    print("markX", d, i_start, i_end)
    for i in range(i_start, i_end):
        j=d-i
        s = mat[i][j]
        explored[i][j] = 1
        if s < max_score - X_drop:
            mat[i][j] = float('-inf')
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
def nw4(A,B, x_thresh=2):
    mat = initialize(A,B, 0.1)
    m = len(A) + 1 # rows
    n = len(B) + 1 # cols
    i_start = 1
    i_end = 1
    max_score = float('-inf')
    # TODO: figure out how to incorporate the base case into the compute and maybe start from 0?
    lower_diag = float('-inf')
    upper_diag = float('inf')
    # this has to be the long side
    for ad in range(2, m): # m is len(A) + 1
        # change to 0 and i_start = 0, i_end = 0 if base case incorporated
        i_end = i_end + 1 # up to but not including
        for i in range(i_start, i_end):
            j=ad-i
            s = score(mat, A, B, i, j)
            mat[i][j] = s
            explored[i][j] = 1
            max_score = max(mat[i][j], max_score)
        if i_start >= i_end:
            break
        markX(mat, ad, i_start, i_end, x_thresh, max_score)
        i_start, d_lo = get_i_start(mat, ad, i_start, i_end)
        i_end, d_hi = get_i_end(mat, ad, i_start, i_end)
        lower_diag = max(lower_diag, d_lo)
        upper_diag = min(upper_diag, d_hi)
        i_start = modify_i_start(i_start, lower_diag, ad)
        i_end = modify_i_end(i_end, upper_diag, ad)
        print(i_start, i_end)
        viz(mat, A, B)
    if i_start >= i_end:
        return mat
    print("----")
    print(ad)
    print("----")
    for ad in range(m, m + n - 1):
        expected_i_start = ad - m + 1
        i_start = max(i_start, expected_i_start)
        expected_i_end = m
        i_end = min(i_end + 1, expected_i_end)
        for i in range(i_start, i_end):
            j=ad-i
            s = score(mat, A, B, i, j)
            mat[i][j] = s
            explored[i][j] = 1
            max_score = max(mat[i][j], max_score)
        if i_start == i_end:
            break
        markX(mat, ad, i_start, i_end, x_thresh, max_score)
        i_start, d_lo = get_i_start(mat, ad, i_start, i_end)
        i_end, d_hi = get_i_end(mat, ad, i_start, i_end)
        lower_diag = max(lower_diag, d_lo)
        upper_diag = min(upper_diag, d_hi)
        i_start = modify_i_start(i_start, lower_diag, ad)
        i_end = modify_i_end(i_end, upper_diag, ad)
        viz(mat, A, B)



    return mat


#A = "ACGGGGAACGGGGAGG"
#B = "ATTCGGAACGGGGAGG"
#B = "ATTCGGAATTCGGA"
mat = nw4(A, B)
print_matrix(mat)