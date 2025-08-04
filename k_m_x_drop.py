def lev_T(a, b):
    m, n = len(a), len(b)
    K = m + n
    # Initialize T with (K+1) × (n+1), filled with large numbers
    T = [[K + 1] * (n + 1) for _ in range(K + 1)]

    # Base cases: j = 0 → T[k][0] = k (i = k); and when i=0 → k=j: T[k][k] = k
    for k in range(K + 1):
        if k <= n:
            T[k][k] = k  # covers the i=0, j=k case
        T[k][0] = k  # covers the j=0 case

    # Fill in DP using your recurrence
    for k in range(1, K + 1):
        for j in range(max((k - 1) // 2, max(1, k - m)), min(n, k) + 1):
            i = k - j
            if i == 0:
                continue  # already set via base T[k][k]
            cost = 0 if a[i - 1] == b[j - 1] else 1
            if cost == 0:
                # match case
                T[k][j] = T[k - 2][j - 1]
            else:
                # mismatch: deletion (k-1,j), insertion (k-1,j-1), substitution (k-2,j-1)
                T[k][j] = 1 + min(T[k - 1][j], T[k - 1][j - 1], T[k - 2][j - 1])

    # Result is S[m][n] = T[m + n][n]
    return T[K][n]


def print_matrix(T, a, b, title="DP Matrix T"):
    """Print the DP matrix in a formatted way"""
    m, n = len(a), len(b)
    K = m + n

    print(f"\n{title}")
    print(f"String a: '{a}' (length {m})")
    print(f"String b: '{b}' (length {n})")
    print(f"K = m + n = {K}")
    print("\nMatrix T[k][j] where k = i + j:")

    # Print column headers
    print("k\\j ", end="")
    for j in range(n + 1):
        print(f"{j:4}", end="")
    print()

    # Print rows
    for k in range(K + 1):
        print(f"{k:2}  ", end="")
        for j in range(n + 1):
            val = T[k][j]
            if val == K + 1:
                print("  - ", end="")
            else:
                print(f"{val:4}", end="")
        print()


def test_lev_T():
    """Test the lev_T function with various examples"""
    print("Testing lev_T function")
    print("=" * 50)

    # Test case 1: Simple example
    a1, b1 = "catttkfopefopwefpoet", "doggggfneifewwfpg"
    result1 = lev_T(a1, b1)
    print(f"\nTest 1: a='{a1}', b='{b1}'")
    print(f"Edit distance: {result1}")

    # Re-run to get matrix for printing
    m, n = len(a1), len(b1)
    K = m + n
    T = [[K + 1] * (n + 1) for _ in range(K + 1)]

    # Base cases
    for k in range(K + 1):
        if k <= n:
            T[k][k] = k
        T[k][0] = k

    # Fill DP
    for k in range(1, K + 1):
        for j in range(max((k - 5) // 2, max(1, k - m)), min(n, k) + 1):
            # i = k - j
            # if i == 0:
            #     continue
            # cost = 0 if a1[i - 1] == b1[j - 1] else 1
            # if cost == 0:
            #     T[k][j] = T[k - 2][j - 1]
            # else:
            #     T[k][j] = 1 + min(T[k - 1][j], T[k - 1][j - 1], T[k - 2][j - 1])
            T[k][j] = 10

    print_matrix(T, a1, b1, "Test 1 Matrix")

    # take matrix T and fill in locations on i, j matrix
    T2 = [[0] * (m + 2) for _ in range(n + 2)]
    for k in range(1, K + 1):
        for j in range(max(1, k - m), min(n, k) + 1):
            i = k - j
            if T[k][j] == 10:
                print(f"({i},{j})", end=" ")
                T2[i][j] = 10

    print("\nMatrix T2 (formatted):")
    for i in range(m + 1):
        for j in range(n + 1):
            if T2[i][j] == 10:
                print(f" 10 ", end="")
            else:
                print("  - ", end="")
        print()

    # Test case 2: Empty strings
    # a2, b2 = "", "abc"
    # result2 = lev_T(a2, b2)
    # print(f"\n\nTest 2: a='{a2}', b='{b2}'")
    # print(f"Edit distance: {result2}")

    # # Test case 3: Identical strings
    # a3, b3 = "hello", "hello"
    # result3 = lev_T(a3, b3)
    # print(f"\nTest 3: a='{a3}', b='{b3}'")
    # print(f"Edit distance: {result3}")

    # # Test case 4: One character difference
    # a4, b4 = "kitten", "sitting"
    # result4 = lev_T(a4, b4)
    # print(f"\nTest 4: a='{a4}', b='{b4}'")
    # print(f"Edit distance: {result4}")

    # print("\n" + "=" * 50)
    # print("All tests completed!")


if __name__ == "__main__":
    test_lev_T()
