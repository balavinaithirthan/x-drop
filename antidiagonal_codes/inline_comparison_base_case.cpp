#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <algorithm>

// Antidiagonal with inline base-case handling
int nw_antidiag_inline(const std::string& A, const std::string& B, int match, int mismatch, int gap,
                       std::vector<std::vector<int>>& dp) {
    int m = A.size(), n = B.size();

    for (int d = 0; d <= m + n; ++d) {
        int i_start = std::max(0, d - n);
        int i_end = std::min(m, d);
        for (int i = i_start; i <= i_end; ++i) {
            int j = d - i;
            if (i == 0 && j == 0) {
                dp[i][j] = 0;
            } else if (i == 0) {
                dp[i][j] = j * gap;
            } else if (j == 0) {
                dp[i][j] = i * gap;
            } else {
                int diag = dp[i - 1][j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch);
                int up   = dp[i - 1][j] + gap;
                int left = dp[i][j - 1] + gap;
                dp[i][j] = std::max({diag, up, left});
            }
        }
    }
    return dp[m][n];
}

// Antidiagonal with preinitialized base-cases
int nw_antidiag_preinit(const std::string& A, const std::string& B, int match, int mismatch, int gap,
                        std::vector<std::vector<int>>& dp) {
    int m = A.size(), n = B.size();
    for (int i = 0; i <= m; ++i) dp[i][0] = i * gap;
    for (int j = 0; j <= n; ++j) dp[0][j] = j * gap;

    for (int d = 2; d <= m + n; ++d) {
        int i_start = std::max(1, d - n);
        int i_end = std::min(m, d - 1);
        for (int i = i_start; i <= i_end; ++i) {
            int j = d - i;
            int diag = dp[i - 1][j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch);
            int up   = dp[i - 1][j] + gap;
            int left = dp[i][j - 1] + gap;
            dp[i][j] = std::max({diag, up, left});
        }
    }
    return dp[m][n];
}

int main() {
    // Generate random DNA sequences of length 1000
    const int len = 10000;
    std::string A(len, 'A'), B(len, 'A');
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 3);
    const char bases[4] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < len; ++i) {
        A[i] = bases[dist(rng)];
        B[i] = bases[dist(rng)];
    }

    int match = 1, mismatch = -1, gap = -2;

    auto run_and_time = [&](auto func, const std::string& name) {
        std::vector<std::vector<int>> dp(len + 1, std::vector<int>(len + 1));

        auto start = std::chrono::high_resolution_clock::now();
        int score = func(A, B, match, mismatch, gap, dp);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        std::cout << name << " score = " << score
                  << "; time = " << elapsed.count() << " s\n";
    };

    run_and_time(nw_antidiag_inline,  "Antidiagonal (inline base-case)");
    run_and_time(nw_antidiag_preinit, "Antidiagonal (preinit base-case)");

    return 0;
}
