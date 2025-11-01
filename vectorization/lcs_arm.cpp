#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <arm_neon.h>

// clang++ -O3 lcs_arm.cpp -o vector_add

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
int LCS_k_m(std::vector<int> &S, const std::string &A, const std::string &B) {

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
            } // LCS doesn't have row and col base cases
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
    // print out matrix
    // std::cout << "LCS_k_m matrix normal->:\n";
    // for (int k = 0; k <= K; ++k)
    // {
    //     for (int j = 0; j <= n; ++j)
    //     {
    //         std::cout << S[idx(k, j)] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    return S[idx(K, n)];
}


int dilated_LCS_NEON(std::vector<int> &S, const std::string &A, const std::string &B)
{
    int n = A.size();
    int m = B.size();

    auto idx = [&](int i, int j) {
        return i * (m + 1) + j;
    };

    for (int i = 0; i <= n; ++i)
    {
        // Handle j < 4 scalar
        for (int j = 0; j < std::min(4, m + 1); ++j)
        {
            if (i == 0 || j == 0) {
                S[idx(i, j)] = 0;
            } else {
                int match = (A[i - 1] == B[j - 1]) ? 1 : 0;
                S[idx(i, j)] = std::max({S[idx(i, j - 1)],
                                         S[idx(i - 1, j)],
                                         S[idx(i - 1, j - 1)] + match});
            }
        }

        // Vectorized j >= 4
        for (int j = 4; j <= m; j += 4)
        {
            // Load S[i][j-4]
            int32x4_t v_sijm4 = vld1q_s32(&S[idx(i, j - 4)]);

            // Load S[i-1][j-4]
            int32x4_t v_si1jm4 = vld1q_s32(&S[idx(i - 1, j - 4)]);

            // Load matches for positions (j-4 ... j-1)
            int32_t match_arr[4];
            for (int lane = 0; lane < 4; lane++) {
                match_arr[lane] = (A[i - 1] == B[(j - 4) + lane]) ? 1 : 0;
            }
            int32x4_t v_match = vld1q_s32(match_arr);

            // Candidate: S[i-1][j-4] + match
            int32x4_t v_cand1 = vaddq_s32(v_si1jm4, v_match);

            // Load the other references S[i-1][j-3], S[i-1][j-2], S[i-1][j-1], S[i-1][j]
            int32_t other_arr[4] = {
                S[idx(i - 1, j - 3)],
                S[idx(i - 1, j - 2)],
                S[idx(i - 1, j - 1)],
                S[idx(i - 1, j)]
            };
            int32x4_t v_other = vld1q_s32(other_arr);

            // Candidate: each of those + match with their own column
            int32_t match_other_arr[4] = {
                other_arr[0] + ((A[i - 1] == B[j - 3]) ? 1 : 0),
                other_arr[1] + ((A[i - 1] == B[j - 2]) ? 1 : 0),
                other_arr[2] + ((A[i - 1] == B[j - 1]) ? 1 : 0),
                other_arr[3] + ((A[i - 1] == B[j])     ? 1 : 0)
            };
            int32x4_t v_other_match = vld1q_s32(match_other_arr);

            // Load S[i-1][j] for last candidate
            int32_t last_arr[4] = {
                S[idx(i - 1, j)], 0, 0, 0
            };
            int32x4_t v_last = vld1q_s32(last_arr);

            // Combine all candidates
            int32x4_t v_max = vmaxq_s32(v_sijm4, v_cand1);
            v_max = vmaxq_s32(v_max, v_other);
            v_max = vmaxq_s32(v_max, v_other_match);
            v_max = vmaxq_s32(v_max, v_last);

            // Store back S[i][j ... j+3]
            vst1q_s32(&S[idx(i, j)], v_max);
        }
    }

    return S[idx(n, m)];
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

int LCS_k_m_neon_tiled(std::vector<int> &S, const std::string &A, const std::string &B)
{
    const int m = (int)A.size();
    const int n = (int)B.size();
    const int K = m + n;

    auto idx = [&](int k, int j) { return k * (n + 1) + j; };

    constexpr int VEC_SIZE = 4;     // NEON 128-bit, 4 x int32 lanes
    constexpr int TILE_J   = 64;    // tune: multiple of 4; try 64/128/256

    for (int k = 0; k <= K; ++k) {
        const int m_start = std::max(0, k - m);
        const int m_end   = std::min(n, k);

        // strip-mine / tile j
        for (int j0 = m_start; j0 <= m_end; j0 += TILE_J) {
            const int j1        = std::min(j0 + TILE_J - 1, m_end);
            const int total     = j1 - j0 + 1;
            const int vec_len   = (total / VEC_SIZE) * VEC_SIZE; // largest multiple of 4 within this tile
            const int j_vec_end = j0 + vec_len;                  // exclusive

            // ===== SIMD body, 4 columns at a time =====
            for (int j = j0; j < j_vec_end; j += VEC_SIZE) {
                // B[j-1..j+2] ascending
                const uint8_t* pB = reinterpret_cast<const uint8_t*>(B.data() + (j - 1));
                const uint8x8_t b8 = vld1_u8(pB);        // low 4 lanes are B[j-1..j+2]
                const uint8x8_t b4 = b8;

                // A[(k-j-1),(k-j-2),(k-j-3),(k-j-4)] descending -> load from (k-j-4) and reverse 64b
                const int a_start = (k - j - 1) - 3;      // k - j - 4
                const uint8_t* pA  = reinterpret_cast<const uint8_t*>(A.data() + a_start);
                const uint8x8_t a8    = vld1_u8(pA);
                const uint8x8_t arev8 = vrev64_u8(a8);    // lanes reordered so low 4 = A[k-j-1],A[k-j-2],A[k-j-3],A[k-j-4]
                const uint8x8_t a4    = arev8;

                // match 0/1 per lane: (a4 == b4)
                const uint8x8_t cmp8   = vceq_u8(a4, b4);         // 0xFF or 0x00
                const uint8x8_t ones8  = vshr_n_u8(cmp8, 7);      // 255>>7 = 1
                const uint16x8_t ones16 = vmovl_u8(ones8);
                const uint32x4_t ones32 = vmovl_u16(vget_low_u16(ones16));
                const int32x4_t  vmatch = vreinterpretq_s32_u32(ones32);

                // prev rows
                const int32x4_t vprev1 = vld1q_s32(&S[idx(k - 1, j    )]); // S[k-1, j..j+3]
                const int32x4_t vprev2 = vld1q_s32(&S[idx(k - 1, j - 1)]); // S[k-1, j-1..j+2]
                const int32x4_t vprev3 = vld1q_s32(&S[idx(k - 2, j - 1)]); // S[k-2, j-1..j+2]

                // prev3 + match, then max across the three
                const int32x4_t vpm    = vaddq_s32(vprev3, vmatch);
                const int32x4_t vmax12 = vmaxq_s32(vprev1, vprev2);
                const int32x4_t vbest  = vmaxq_s32(vmax12, vpm);

                vst1q_s32(&S[idx(k, j)], vbest);
            }
        }
    }

    return S[idx(K, n)];
}


// Diagonal-style LCS with (k, m) indexing, where k = i + j and m = j
int LCS_k_m_neon(std::vector<int> &S, const std::string &A, const std::string &B)
{
    int m = A.size();
    int n = B.size();
    int K = m + n;

    auto idx = [&](int k, int j)
    { return k * (n + 1) + j; };
    constexpr int VEC_SIZE = 4; // Vector size for NEON
    for (int k = 0; k <= K; ++k)
    {

        int m_start = std::max(0, k - m);
        int m_end = std::min(n, k);
        int total = m_end - m_start + 1;
        int vec_len = (total/ VEC_SIZE) * 4;
        
        // j runs in steps of 4
        for (int j = m_start; j < m_start + vec_len; j += 4)
        {
            // j vector: [j, j+1, j+2, j+3]
            // int32x4_t vj = vaddq_s32(vdupq_n_s32(j), (int32x4_t){0, 1, 2, 3}); // create a vector of 4 ints, [j, j+1, j+2, j+3]
            // int32x4_t vk = vdupq_n_s32(k); // create a vector of 4 ints, [k, k, k, k]
            // int32x4_t vkj = vsubq_s32(vk, vj); // k - j per lane, [k-j, k-j-1, k-j +-2, k-j + -3]

            // Base-case mask: (j == 0) OR (k - j == 0)
            // uint32x4_t mask_j0 = vceqq_s32(vj, vdupq_n_s32(0)); // compare vj with 0, [j==0, j+1==0, j+2==0, j+3==0]
            // uint32x4_t mask_kj0 = vceqq_s32(vkj, vdupq_n_s32(0)); // compare vkj with 0, [k-j==0, k-j-1==0, k-j-2==0, k-j-3==0]
            // uint32x4_t mask_base = vorrq_u32(mask_j0, mask_kj0); // OR the two masks, [j==0 || k-j==0, j+1==0 || k-j-1==0, j+2==0 || k-j-2==0, j+3==0 || k-j-3==0]

            // Load B[j-1..j+2] as bytes (contiguous increasing)
            const uint8_t *pB = reinterpret_cast<const uint8_t *>(B.data() + (j - 1));
            uint8x8_t b8 = vld1_u8(pB); // we only need the low 4 lanes
            uint8x8_t b8lo = b8;

            // Load A[(k-j-1),(k-j-2),(k-j-3),(k-j-4)] as bytes.
            // These are contiguous *decreasing*, so load from the lowest address and reverse.
            int a_start = (k - j - 1) - 3; // equals k - j - 4
            const uint8_t *pA = reinterpret_cast<const uint8_t *>(A.data() + a_start);
            uint8x8_t a8 = vld1_u8(pA); // load 8 to be safe
            // Reverse the first 4 bytes into ascending-lane order
            // Do a 64-bit reverse then we’ll use the low 4 lanes.
            uint8x8_t a8rev = vrev64_u8(a8); // reverse bytes within 64-bit halves
            // Now a8rev low 4 lanes correspond to A[k-j-1], A[k-j-2], A[k-j-3], A[k-j-4]
            uint8x16_t a16 = vcombine_u8(a8, vdup_n_u8(0)); // upper half zeroed
            uint8x16_t b16 = vcombine_u8(b8, vdup_n_u8(0));
            uint8x16_t cmp16 = vceqq_u8(a16, b16); // ✔ OK now

            uint16x4_t cmp16lo = vget_low_u16(cmp16); // lower 4 lanes
            uint32x4_t cmp32u = vmovl_u16(cmp16lo);
            // cmp32u lanes are 0 or 0x000000FF; turn 0x..FF into 1 by compare-with-zero
            uint32x4_t ones_mask = vcgtq_u32(cmp32u, vdupq_n_u32(0));
            int32x4_t vmatch = vandq_s32(vreinterpretq_s32_u32(ones_mask), vdupq_n_s32(1));

            // Load prevs: each is 4 contiguous ints
            int32x4_t vprev1 = vld1q_s32(&S[idx(k - 1, j)]);     // S[k-1, j..j+3]
            int32x4_t vprev2 = vld1q_s32(&S[idx(k - 1, j - 1)]); // S[k-1, j-1..j+2]
            int32x4_t vprev3 = vld1q_s32(&S[idx(k - 2, j - 1)]); // S[k-2, j-1..j+2]

            // prev3 + match
            int32x4_t vpm = vaddq_s32(vprev3, vmatch);

            // max(prev1, prev2, prev3+match)
            int32x4_t vmax12 = vmaxq_s32(vprev1, vprev2);
            int32x4_t vbest = vmaxq_s32(vmax12, vpm); // 

            // If base mask, write 0; else write vbest
            // int32x4_t vbest = vdupq_n_s32(1); // create a vector of 4 ints, [1, 1, 1, 1]
            // int32x4_t vzeros = vdupq_n_s32(0); // create a vector of 4 ints, [0, 0, 0, 0]
            // int32x4_t vout = vbslq_s32(mask_base, vzeros, vbest); // if mask_base is true, write 0, else write vbest

            vst1q_s32(&S[idx(k, j)], vbest);
        }
        // for (int j = m_start + vec_len; j <= m_end; ++j)
        // {
        //     int i = k - j;
        //     if (i == 0 || j == 0)
        //     {
        //         S[idx(k, j)] = 0;
        //     } // LCS doesn't have row and col base cases
        //     else
        //     {
        //         int prev1 = S[idx(k - 1, j)];
        //         int prev2 = S[idx(k - 1, j - 1)];
        //         int prev3 = S[idx(k - 2, j - 1)];
        //         int match = (A[i - 1] == B[j - 1]) ? 1 : 0;
        //         S[idx(k, j)] = std::max({prev1, prev2, prev3 + match});
        //     }
        // }
    }
    // print out matrix
    // std::cout << "LCS_k_m matrix normal->:\n";
    // for (int k = 0; k <= K; ++k)
    // {
    //     for (int j = 0; j <= n; ++j)
    //     {
    //         std::cout << S[idx(k, j)] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    return S[idx(K, n)];
}

void neon_vs_regular(std::string A, std::string B)
{
    double total_elapsed_neon = 0;
    double total_elapsed_regular = 0;
    double total_elapsed_dilated_neon = 0;
    for (int i = 0; i < 10; ++i)
    {
        int n = A.size();
        int m = B.size();
        int K = n + m;
        std::vector<int> S_km_neon((K + 1) * (m + 1), 0);
        std::vector<int> S((n + 1) * (m + 1), 0);
        std::vector<int> S_dilated_neon((n + 1) * (m + 1), 0);
        auto start = std::chrono::high_resolution_clock::now();
        int result5 = LCS_k_m_neon(S_km_neon, A, B);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        total_elapsed_neon += elapsed.count();
        start = std::chrono::high_resolution_clock::now();
        int result6 = LCS(A, B, S);
        end = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        total_elapsed_regular += elapsed.count();
        start = std::chrono::high_resolution_clock::now();
        int result7 = dilated_LCS_NEON(S_dilated_neon, A, B);
        end = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        total_elapsed_dilated_neon += elapsed.count();
    }
    std::cout << "Total elapsed time for NEON: " << total_elapsed_neon << " milliseconds\n";
    std::cout << "Total elapsed time for regular: " << total_elapsed_regular << " milliseconds\n";
    std::cout << "Total elapsed time for dilated NEON: " << total_elapsed_dilated_neon << " milliseconds\n";
    std::cout << "Ratio: " << total_elapsed_regular/total_elapsed_neon  << std::endl;
    std::cout << "Ratio: " << total_elapsed_regular / total_elapsed_dilated_neon << std::endl;
}

int main()
{
    std::string A(4000, 'A'); // 1000 As
    std::string B(4000, 'A'); // 1000 As
    // std::string A = "BACACC";
    // std::string B = "ABCBAC";

    int n = A.size();
    int m = B.size();
    int K = n + m;

    // Standard LCS
    std::vector<int> dp((n + 1) * (m + 1), 0);
    int result1 = LCS(A, B, dp);
    std::cout << "Classic (i,j) LCS result: " << result1 << std::endl;

    std::cout << "______________________: " << std::endl;
    // Diagonal LCS using (k, m) coordinates
    std::vector<int> S_km((K + 1) * (m + 1), 0);
    int result3 = LCS_k_m(S_km, A, B);
    std::cout << "Diagonal (k,m) LCS result: " << result3 << std::endl;
    std::cout << "______________________: " << std::endl;

    // Dilated LCS
    std::vector<int> S_dilated((n + 1) * (m + 1), 0);
    int result4 = dilated_LCS(S_dilated, A, B);
    std::cout << "Dilated LCS result: " << result4 << std::endl;
    std::cout << "______________________: " << std::endl;

    // Dilated LCS with NEON

    // Diagonal LCS using (k, m) coordinates with NEON
    std::vector<int> S_km_neon((K + 1) * (m + 1), 0);
    int result5 = LCS_k_m_neon(S_km_neon, A, B);
    std::cout << "Diagonal (k,m) LCS result with NEON: " << result5 << std::endl;

    std::cout << "______________________: " << std::endl;
    std::vector<int> S_km_neon_tiled((K + 1) * (m + 1), 0);
    int result6 = LCS_k_m_neon_tiled(S_km_neon_tiled, A, B);
    std::cout << "Diagonal (k,m) LCS result with NEON TILED: " << result6 << std::endl;
    std::cout << "______________________: " << std::endl;

    // repeat neon and regular LCS 10 times and print the average time
    neon_vs_regular(A, B);

    return 0;
}
