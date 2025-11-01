import numpy as np

# Input strings
A = "kitten"
B = "sitting"
m = len(A)
n = len(B)

# Initialize the score-indexed DP matrix
T = np.full((m + n + 1, n + 1), np.inf)
visited = np.zeros((m + n + 1, n + 1), dtype=bool)

# Base case: when j = 0, i = k, so cost is i
for i in range(m + 1):
    T[i][0] = i


# Fill the DP matrix using the transformed recurrence
for k in range(0, m + n + 1):
    j_min = max(1, k - m)
    j_max = min(n, k)
    for j in range(j_min, j_max + 1):
        i = k - j
        if not (0 <= i < m):
            continue

        insertion = T[k - 1][j - 1] + 1
        deletion = T[k - 1][j] + 1
        substitution = T[k - 2][j - 1] + (A[i] != B[j - 1])
        visited[k][j] = True

        T[k][j] = min(insertion, deletion, substitution)


# Format the matrix for printing
def format_matrix_aligned(matrix, width=4):
    return [
        [
            "{:<{w}}".format(int(cell) if np.isfinite(cell) else "inf", w=width)
            for cell in row
        ]
        for row in matrix
    ]


# Print the matrix
print("Score-Indexed DP Matrix (T[k][j] with k = i + j):")
print(" " * 6 + "  ".join(f"j={j}" for j in range(n + 1)))
for k, row in enumerate(format_matrix_aligned(T)):
    print(f"k={k:<2}  " + "  ".join(row))

print(" " * 6 + "  ".join(f"j={j}" for j in range(n + 1)))
for k, row in enumerate(format_matrix_aligned(visited)):
    print(f"k={k:<2}  " + "  ".join(row))
