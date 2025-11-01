def lev_standard(a, b):
    m, n = len(a), len(b)
    S = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        S[i][0] = i
    for j in range(n + 1):
        S[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            S[i][j] = min(S[i - 1][j] + 1, S[i][j - 1] + 1, S[i - 1][j - 1] + cost)
    return S[m][n]


def lev_wavefront(a, b):
    m, n = len(a), len(b)
    S = [[0] * (n + 1) for _ in range(m + 1)]
    for k in range(m + n + 1):
        for j in range(max(0, k - m), min(n, k) + 1):
            i = k - j
            if i == 0 or j == 0:
                S[i][j] = i + j
            else:
                cost = 0 if a[i - 1] == b[j - 1] else 1
                S[i][j] = min(S[i - 1][j] + 1, S[i][j - 1] + 1, S[i - 1][j - 1] + cost)
    return S[m][n]


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
        for j in range(max(1, k - m), min(n, k) + 1):
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


# 🔍 Compare outputs and timing
if __name__ == "__main__":
    tests = [("kitten", "sitting"), ("flaw", "lawn"), ("gumbo", "gambol")]
    import time

    for a, b in tests:
        results = {}
        for name, fn in [
            ("standard", lev_standard),
            ("wavefront", lev_wavefront),
            ("two_rows", lev_T),
        ]:
            start = time.time()
            d = fn(a, b)
            elapsed = time.time() - start
            results[name] = (d, elapsed)
        print(f"\nTest: {a!r} → {b!r}")
        for name, (d, t) in results.items():
            print(f"  {name:10}: dist={d}, time={t*1e6:.1f} µs")
