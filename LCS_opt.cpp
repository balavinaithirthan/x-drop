#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <immintrin.h>

// Standard LCS with (i, j) indexing
int LCS(const std::string &A, const std::string &B, std::vector<int> &dp)
{
    int n = A.size();
    int m = B.size();

    auto idx = [&](int i, int j)
    { return i * (m + 1) + j; };

    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= m; ++j)
        {
            if (i == 0 || j == 0)
            {
                dp[idx(i, j)] = 0;
            }
            else if (A[i - 1] == B[j - 1])
            {
                dp[idx(i, j)] = dp[idx(i - 1, j - 1)] + 1;
            }
            else
            {
                dp[idx(i, j)] = std::max(dp[idx(i - 1, j)], dp[idx(i, j - 1)]);
            }
        }
    }

    return dp[idx(n, m)];
}

// Diagonal-style LCS with (k, m) indexing, where k = i + j and m = j
int LCS_k_m(std::vector<int> &S, const std::string &A, const std::string &B)
{
    int m = A.size();
    int n = B.size();
    int K = m + n;

    auto idx = [&](int k, int j)
    { return k * (n + 1) + j; };

    for (int k = 0; k <= K; ++k)
    {
        int m_start = std::max(0, k - m);
        int m_end = std::min(n, k);

        for (int j = m_start; j <= m_end; ++j)
        {
            int i = k - j;
            if (i == 0 || j == 0)
            {
                S[idx(k, j)] = 0;
            }
            else
            {
                int prev1 = S[idx(k - 1, j)];
                int prev2 = S[idx(k - 1, j - 1)];
                int prev3 = S[idx(k - 2, j - 1)];
                int match = (A[i - 1] == B[j - 1]) ? 1 : 0;
                S[idx(k, j)] = std::max({prev1, prev2, prev3 + match});
            }
        }
    }

    return S[idx(K, n)];
}

int dilated_LCS(std::vector<int> &S, const std::string &A, const std::string &B)
{
    int n = A.size();
    int m = B.size();

    auto idx = [&](int i, int j)
    {
        return i * (m + 1) + j;
    };

    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= m; ++j)
        {
            if (i == 0 || j == 0)
            {
                S[idx(i, j)] = 0;
            }
            else if (j < 4)
            {
                int match = (A[i - 1] == B[j - 1]) ? 1 : 0;
                S[idx(i, j)] = std::max({S[idx(i, j - 1)],
                                         S[idx(i - 1, j)],
                                         S[idx(i - 1, j - 1)] + match});
            }
            else
            {
                S[idx(i, j)] = std::max({S[idx(i, j - 4)],
                                         S[idx(i - 1, j - 4)] + (A[i - 1] == B[j - 4]),
                                         S[idx(i - 1, j - 3)],
                                         S[idx(i - 1, j - 3)] + (A[i - 1] == B[j - 3]),
                                         S[idx(i - 1, j - 2)],
                                         S[idx(i - 1, j - 2)] + (A[i - 1] == B[j - 2]),
                                         S[idx(i - 1, j - 1)],
                                         S[idx(i - 1, j - 1)] + (A[i - 1] == B[j - 1]),
                                         S[idx(i - 1, j)]});
            }
        }
    }

    return S[idx(n, m)];
}

// Computes LCS via anti-diagonals (k = i + j), vectorizing the inner j-loop by 8.
int LCS_k_m_vec(std::vector<int> &S,
                const std::string &A,
                const std::string &B)
{
    const int m = A.size();
    const int n = B.size();
    const int K = m + n;
    S.assign((K + 1) * (n + 1), 0);

    // Pre-convert B into 32-bit ints for easy gathering
    std::vector<int> B_int(n + 1);
    for (int j = 1; j <= n; ++j)
    {
        B_int[j] = B[j - 1];
    }

    auto idx = [&](int k, int j)
    { return k * (n + 1) + j; };
    const int W = 8; // vector width
    const __m256i one = _mm256_set1_epi32(1);

    for (int k = 0; k <= K; ++k)
    {
        int j0 = std::max(0, k - m);
        int j1 = std::min(n, k);

        int j = j0;
        while (j <= j1)
        {
            // check if we can do a full 8-wide vector step
            if (j + W - 1 <= j1 && j > 0 && (k > 0) && (k > 1))
            {
                int i = k - j; // for match, but we'll build a broadcast per lane below

                // load prev1 = S[k-1][j .. j+7]
                __m256i prev1 = _mm256_loadu_si256(
                    (const __m256i *)&S[idx(k - 1, j)]);

                // load prev2 = S[k-1][j-1 .. j+6]
                __m256i prev2 = _mm256_loadu_si256(
                    (const __m256i *)&S[idx(k - 1, j - 1)]);

                // load prev3 = S[k-2][j-1 .. j+6]
                __m256i prev3 = _mm256_loadu_si256(
                    (const __m256i *)&S[idx(k - 2, j - 1)]);

                // gather B[j-1 .. j+6] into a vector
                __m256i idxs = _mm256_setr_epi32(
                    j - 1, j, j + 1, j + 2,
                    j + 3, j + 4, j + 5, j + 6);
                __m256i Bv = _mm256_i32gather_epi32(
                    B_int.data(), idxs, 4);

                // broadcast A[i-1]
                __m256i Av = _mm256_set1_epi32(A[i - 1]);

                // match = (A[i-1] == B[j-1..j+6]) ? 1 : 0
                __m256i eq = _mm256_cmpeq_epi32(Av, Bv);
                __m256i match = _mm256_and_si256(eq, one);

                // prev3 + match
                __m256i prev3m = _mm256_add_epi32(prev3, match);

                // tmp = max(prev2, prev3+match)
                __m256i tmp = _mm256_max_epi32(prev2, prev3m);

                // result = max(prev1, tmp)
                __m256i res = _mm256_max_epi32(prev1, tmp);

                // store back into S[k][j .. j+7]
                _mm256_storeu_si256(
                    (__m256i *)&S[idx(k, j)], res);

                j += W;
            }
            else
            {
                // Scalar fallback for head, tail, or boundary conditions
                int i = k - j;
                if (i == 0 || j == 0)
                {
                    S[idx(k, j)] = 0;
                }
                else
                {
                    int p1 = S[idx(k - 1, j)];
                    int p2 = S[idx(k - 1, j - 1)];
                    int p3 = S[idx(k - 2, j - 1)];
                    int match = (A[i - 1] == B[j - 1]) ? 1 : 0;
                    S[idx(k, j)] = std::max({p1, p2, p3 + match});
                }
                ++j;
            }
        }
    }

    return S[idx(K, n)];
}

int main()
{
    std::string A = "AGGTABCGGTTAGCGATC";
    std::string B = "GXTXAYBGGCATGCTAGC";

    int n = A.size();
    int m = B.size();
    int K = n + m;

    // Standard LCS
    std::vector<int> dp((n + 1) * (m + 1), 0);
    int result1 = LCS(A, B, dp);
    std::cout << "Classic (i,j) LCS result: " << result1 << std::endl;

    // Diagonal LCS using (k, m) coordinates
    std::vector<int> S_km((K + 1) * (m + 1), 0);
    int result2 = LCS_k_m(S_km, A, B);
    std::cout << "Diagonal (k,m) LCS result: " << result2 << std::endl;

    // vectorized LCS using (k,m) coordinates
    std::vector<int> S_km((K + 1) * (m + 1), 0);
    int result2 = LCS_k_m_vec(S_km, A, B);
    std::cout << "Diagonal (k,m) LCS result: " << result2 << std::endl;

    // Dilated LCS
    std::vector<int> S_dilated((n + 1) * (m + 1), 0);
    int result3 = dilated_LCS(S_dilated, A, B);
    std::cout << "Dilated LCS result: " << result3 << std::endl;

    return 0;
}
