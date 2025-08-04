def delta(x, y, z):
    gap = "-"
    return int(x != y or y != z or x != z)


def triple_edit_distance(A, B, C):
    m, n, p = len(A), len(B), len(C)
    D = [[[float("inf")] * (p + 1) for _ in range(n + 1)] for _ in range(m + 1)]
    D[0][0][0] = 0

    for i in range(m + 1):
        for j in range(n + 1):
            for k in range(p + 1):
                if i > 0 and j > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k],
                        D[i - 1][j - 1][k - 1] + delta(A[i - 1], B[j - 1], C[k - 1]),
                    )
                if i > 0 and j > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j - 1][k] + delta(A[i - 1], B[j - 1], "-")
                    )
                if i > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j][k - 1] + delta(A[i - 1], "-", C[k - 1])
                    )
                if j > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j - 1][k - 1] + delta("-", B[j - 1], C[k - 1])
                    )
                if i > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j][k] + delta(A[i - 1], "-", "-")
                    )
                if j > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j - 1][k] + delta("-", B[j - 1], "-")
                    )
                if k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j][k - 1] + delta("-", "-", C[k - 1])
                    )
    # print all elements of D[m][n][p] in a formatted way
    print("Triple Edit Distance Matrix (D[i][j][k]):")
    for i in range(m + 1):
        for j in range(n + 1):
            for k in range(p + 1):
                print(f"{D[i][j][k]}", end=" ")
            print()
    return D[m][n][p]


def antidiag_edit_distance(A, B, C):
    m, n, p = len(A), len(B), len(C)
    D = [[[float("inf")] * (p + 1) for _ in range(n + 1)] for _ in range(m + 1)]
    D[0][0][0] = 0

    for k in range(m + n + p + 1):
        for j in range(n + 1):
            for k in range(p + 1):
                if i > 0 and j > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k],
                        D[i - 1][j - 1][k - 1] + delta(A[i - 1], B[j - 1], C[k - 1]),
                    )
                if i > 0 and j > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j - 1][k] + delta(A[i - 1], B[j - 1], "-")
                    )
                if i > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j][k - 1] + delta(A[i - 1], "-", C[k - 1])
                    )
                if j > 0 and k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j - 1][k - 1] + delta("-", B[j - 1], C[k - 1])
                    )
                if i > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i - 1][j][k] + delta(A[i - 1], "-", "-")
                    )
                if j > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j - 1][k] + delta("-", B[j - 1], "-")
                    )
                if k > 0:
                    D[i][j][k] = min(
                        D[i][j][k], D[i][j][k - 1] + delta("-", "-", C[k - 1])
                    )
    # print all elements of D[m][n][p] in a formatted way
    print("Triple Edit Distance Matrix (D[i][j][k]):")
    for i in range(m + 1):
        for j in range(n + 1):
            for k in range(p + 1):
                print(f"{D[i][j][k]}", end=" ")
            print()
    return D[m][n][p]


triple_edit_distance(A="kitten", B="sitting", C="kittening")
antidiag_edit_distance(A="kitten", B="sitting", C="kittening")
