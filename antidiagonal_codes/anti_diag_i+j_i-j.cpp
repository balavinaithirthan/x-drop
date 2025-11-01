#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>

int main()
{
    // Input strings
    std::string q = "kitten";
    std::string t = "sitting";

    int q_len = q.size();
    int t_len = t.size();
    int k_max = q_len + t_len;

    // Matrix E[k][m] → offset m by m_offset so we can index negative m
    int m_range = 2 * k_max + 1;
    int m_offset = k_max; // m ∈ [-k_max, k_max] → shift by k_max

    // E[m][k] with m ∈ [-k_max, k_max], k ∈ [0, k_max]
    std::vector<std::vector<int>> E(m_range, std::vector<int>(k_max + 1, 100));

    // Delta function: 0 if equal, 1 otherwise
    auto delta = [](char a, char b) -> int
    {
        return (a == b) ? 0 : 1;
    };

    for (int k = 0; k <= k_max; ++k)
    {
        for (int m = -k; m <= k; m += 2)
        {
            int m_idx = m + m_offset;

            // Base cases
            if (m == -k || m == k)
            {
                E[m_idx][k] = k;
                continue;
            }

            // Check parity — both i and j must be integers
            if ((k - m) % 2 != 0)
                continue;

            int i = (k - m) / 2;
            int j = (k + m) / 2;

            // Only compute if i and j are valid

            int insertion = E[m + 1 + m_offset][k - 1] + 1;                // (i+1, j)
            int deletion = E[m - 1 + m_offset][k - 1] + 1;                 // (i, j+1)
            int substitution = E[m + m_offset][k - 2] + delta(q[i], t[j]); // (i+1, j+1)

            E[m_idx][k] = std::min({insertion, deletion, substitution});
        }
    }

    // Final result: edit distance between q and t
    int final_k = q_len + t_len;
    int final_m = t_len - q_len;
    int m_idx = final_m + m_offset;

    // Search backwards to find minimal valid E(m, k)
    int result = INT_MAX;
    for (int k = 0; k <= k_max; ++k)
    {
        int m = t_len - q_len;
        int i = (k - m) / 2;
        int j = (k + m) / 2;
        if ((k - m) % 2 == 0 && i == q_len && j == t_len)
        {
            result = std::min(result, E[m + m_offset][k]);
        }
    }
    // print out the matrix
    for (int k = 0; k <= k_max; ++k)
    {
        for (int m = -k_max; m <= k_max; m += 2)
        {
            int m_idx = m + m_offset;
            std::cout << E[m_idx][k] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Edit distance between \"" << q << "\" and \"" << t << "\" is: " << result << std::endl;

    return 0;
}
