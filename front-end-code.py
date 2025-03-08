# GLOBALS
def Matrix():
    pass

matrix = ""

@recuma
def x_drop_nw(seed):
    inputs = matrix, seed
    A = Matrix(length = matrix.X - seed.X, width = matrix.Y - seed.Y)
    Substitution_matrix = {
    ('A', 'A'): 1, ('A', 'T'): -1, ('A', 'C'): -1, ('A', 'G'): -1,
    ('T', 'A'): -1, ('T', 'T'): 1, ('T', 'C'): -1, ('T', 'G'): -1,
    ('C', 'A'): -1, ('C', 'T'): -1, ('C', 'C'): 1, ('C', 'G'): -1,
    ('G', 'A'): -1, ('G', 'T'): -1, ('G', 'C'): -1, ('G', 'G'): 1,
}
    recurrence = {
        S[i][0] = i * G
        S[0][i] = i * G
        S[i][j] = max(S[i-1][j-1] + Substitution_matrix[i][j], S[i-1][j] + G, S[i][j-1] + G)
        max_score = max(max_score, S[i][j])
        max_ij = argmax(max_score, S[i][j])
    }

    mapping = Mapping("antiDiagonal", ad = i + j)
    mapping = Mapping("diagonal", d = (i + j) // 2)
    frontier = Frontier(mapping = ad, index = d, size = 2)
    frontier.reccurence = {
        frontier[d] =  if max_score > frontier[d] - 15 -> frontier[d] = -inf :
        d_lo[t] =  argmin(d, frontier[d] != -inf)
        d_hi[t] = argmin(d, frontier[d] != -inf)
        d_lo[t] = max(d_lo[t], d_lo[t-1])
        d_hi[t] = min(d_hi[t], d_hi[t-1])
    }
    frontier.stop(d_lo > d_hi)
    frontier.optimize(

    )


    outputs = S[-1][-1], max_score, max_ij


def ungapped_filter(seed):
    inputs = 
    mapping = Mapping("antiDiagonal", ad = i + j)
    mapping = Mapping("diagonal", d = (i + j) // 2)
    A = Matrix(length = matrix.X - seed.X, width = matrix.Y - seed.Y)
    recurrence = {
        A[d][ad] = A[d][ad-1] + Substitution_matrix[i][j]
    }
    frontier = Frontier(mapping = ad, index = d, size = 2)
    frontier.reccurence = {
        d_lo = -1
        d_hi = 1
    }
    frontier.stop(Substitution_matrix[i][j] == -1)
    frontier.optimize(

    )

    outputs = 