#include <iostream>
#include <chrono>
#include <string>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>
using std::string;

int a[5001][5001];
int b[8];

class Solution
{
public:
    virtual string name() = 0;
    virtual int longestCommonSubsequence(string text1, string text2) = 0;
    virtual ~Solution() = default;
};

// 16ms
class SolutionNaive : public Solution
{
public:
    string name()
    {
        return "naive solution";
    }
    int longestCommonSubsequence(string text1, string text2)
    {
        int b = (text1[0] == text2[0]);
        a[0][0] = b; // base case L[0][0] = text1[0] == text2[0]
        for (int i = 1; i < text1.length(); i++)
        {
            a[i][0] = std::max(a[i - 1][0], int(text1[i] == text2[0]));
            // L[i][0] = max(L[i-1][0], text1)
        }
        for (int i = 1; i < text2.length(); i++)
        {
            a[0][i] = std::max(a[0][i - 1], int(text1[0] == text2[i]));
            //
        }
        for (int i = 1; i < text1.length(); i++)
        {
            for (int j = 1; j < text2.length(); j++)
            {
                a[i][j] = std::max(std::max(a[i - 1][j], a[i][j - 1]), a[i - 1][j - 1] + (text1[i] == text2[j]));
            }
        }
        return a[text1.length() - 1][text2.length() - 1];
    }
};

// 11ms
class SolutionCompress : public Solution
{
public:
    string name()
    {
        return "compressed solution";
    }
    // only use first and second row
    __attribute__((target("avx2"))) int longestCommonSubsequence(string text1, string text2)
    {
        int b = (text1[0] == text2[0]);
        a[0][0] = b;
        // init base case for each column, row 0
        for (int i = 1; i < text2.length(); i++)
            a[0][i] = std::max(a[0][i - 1], int(text1[0] == text2[i]));
        for (int i = 1; i < text1.length(); i++)
        {
            // move base case inside loop
            a[i & 1][0] = std::max(a[!(i & 1)][0], int(text1[i] == text2[0]));
            for (int j = 1; j < text2.length(); j++)
            {
                a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), a[!(i & 1)][j - 1] + (text1[i] == text2[j]));
            }
        }
        return a[(text1.length() - 1) & 1][text2.length() - 1];
    }
};

// 7ms

class SolutionBy4 : public Solution
{
public:
    string name()
    {
        return "solution vectorized by 4";
    }

    int longestCommonSubsequence(string text1, string text2)
    {
        a[0][0] = int(text1[0] == text2[0]);
        for (int i = 1; i < text2.length(); i++)
            a[0][i] = std::max(a[0][i - 1], int(text1[0] == text2[i]));
        for (int i = 1; i < text1.length(); i++)
        {
            if (text2.length() <= 3)
            {
                a[i & 1][0] = std::max(a[!(i & 1)][0], int(text1[i] == text2[0]));
                for (int j = 1; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
                }
                continue;
            }

            b[0] = int(text1[i] == text2[0]);
            b[1] = a[!(i & 1)][0] + int(text1[i] == text2[1]);
            b[2] = a[!(i & 1)][1] + int(text1[i] == text2[2]);
            b[3] = a[!(i & 1)][2] + int(text1[i] == text2[3]);

            b[1] = std::max(b[1], b[0]);
            b[3] = std::max(b[3], b[2]);

            b[2] = std::max(b[2], b[1]);
            b[3] = std::max(b[3], b[1]);

            a[i & 1][0] = std::max(a[!(i & 1)][0], b[0]);
            a[i & 1][1] = std::max(a[!(i & 1)][1], b[1]);
            a[i & 1][2] = std::max(a[!(i & 1)][2], b[2]);
            a[i & 1][3] = std::max(a[!(i & 1)][3], b[3]);

            int j = 4;
            while (j + 4 <= text2.length())
            {
                int k = b[3];
                b[0] = std::max(k, a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
                b[1] = std::max(k, a[!(i & 1)][j + 0] + int(text1[i] == text2[j + 1]));
                b[2] = std::max(k, a[!(i & 1)][j + 1] + int(text1[i] == text2[j + 2]));
                b[3] = std::max(k, a[!(i & 1)][j + 2] + int(text1[i] == text2[j + 3]));

                b[1] = std::max(b[1], b[0]);
                b[3] = std::max(b[3], b[2]);

                b[2] = std::max(b[2], b[1]);
                b[3] = std::max(b[3], b[1]);

                a[i & 1][j + 0] = std::max(a[!(i & 1)][j + 0], b[0]);
                a[i & 1][j + 1] = std::max(a[!(i & 1)][j + 1], b[1]);
                a[i & 1][j + 2] = std::max(a[!(i & 1)][j + 2], b[2]);
                a[i & 1][j + 3] = std::max(a[!(i & 1)][j + 3], b[3]);

                j += 4;
            }

            for (; j < text2.length(); j++)
            {
                a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
            }
        }
        return a[(text1.length() - 1) & 1][text2.length() - 1];
    }
};

// 11ms
class SolutionBy8 : public Solution
{
public:
    string name()
    {
        return "solution vectorized by 8";
    }

    int longestCommonSubsequence(string text1, string text2)
    {
        a[0][0] = int(text1[0] == text2[0]);
        for (int i = 1; i < text2.length(); i++)
            a[0][i] = std::max(a[0][i - 1], int(text1[0] == text2[i]));

        if (text2.length() <= 8)
        {
            for (int i = 1; i < text1.length(); i++)
            {
                a[i & 1][0] = std::max(a[!(i & 1)][0], int(text1[i] == text2[0]));
                for (int j = 1; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
                }
            }
        }
        else
        {
            for (int i = 1; i < text1.length(); i++)
            {

                b[0] = int(text1[i] == text2[0]);
                b[1] = a[!(i & 1)][0] + int(text1[i] == text2[1]);
                b[2] = a[!(i & 1)][1] + int(text1[i] == text2[2]);
                b[3] = a[!(i & 1)][2] + int(text1[i] == text2[3]);
                b[4] = a[!(i & 1)][3] + int(text1[i] == text2[4]);
                b[5] = a[!(i & 1)][4] + int(text1[i] == text2[5]);
                b[6] = a[!(i & 1)][5] + int(text1[i] == text2[6]);
                b[7] = a[!(i & 1)][6] + int(text1[i] == text2[7]);

                b[1] = std::max(b[1], b[0]);
                b[3] = std::max(b[3], b[2]);
                b[5] = std::max(b[5], b[4]);
                b[7] = std::max(b[7], b[6]);

                b[2] = std::max(b[2], b[1]);
                b[3] = std::max(b[3], b[1]);
                b[6] = std::max(b[6], b[5]);
                b[7] = std::max(b[7], b[5]);

                b[4] = std::max(b[4], b[3]);
                b[5] = std::max(b[5], b[3]);
                b[6] = std::max(b[6], b[3]);
                b[7] = std::max(b[7], b[3]);

                a[i & 1][0] = std::max(a[!(i & 1)][0], b[0]);
                a[i & 1][1] = std::max(a[!(i & 1)][1], b[1]);
                a[i & 1][2] = std::max(a[!(i & 1)][2], b[2]);
                a[i & 1][3] = std::max(a[!(i & 1)][3], b[3]);
                a[i & 1][4] = std::max(a[!(i & 1)][4], b[4]);
                a[i & 1][5] = std::max(a[!(i & 1)][5], b[5]);
                a[i & 1][6] = std::max(a[!(i & 1)][6], b[6]);
                a[i & 1][7] = std::max(a[!(i & 1)][7], b[7]);

                int j = 8;
                while (j + 8 <= text2.length())
                {
                    int k = b[7];
                    b[0] = std::max(k, a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
                    b[1] = std::max(k, a[!(i & 1)][j + 0] + int(text1[i] == text2[j + 1]));
                    b[2] = std::max(k, a[!(i & 1)][j + 1] + int(text1[i] == text2[j + 2]));
                    b[3] = std::max(k, a[!(i & 1)][j + 2] + int(text1[i] == text2[j + 3]));
                    b[4] = std::max(k, a[!(i & 1)][j + 3] + int(text1[i] == text2[j + 4]));
                    b[5] = std::max(k, a[!(i & 1)][j + 4] + int(text1[i] == text2[j + 5]));
                    b[6] = std::max(k, a[!(i & 1)][j + 5] + int(text1[i] == text2[j + 6]));
                    b[7] = std::max(k, a[!(i & 1)][j + 6] + int(text1[i] == text2[j + 7]));

                    b[1] = std::max(b[1], b[0]);
                    b[3] = std::max(b[3], b[2]);
                    b[5] = std::max(b[5], b[4]);
                    b[7] = std::max(b[7], b[6]);

                    b[2] = std::max(b[2], b[1]);
                    b[3] = std::max(b[3], b[1]);
                    b[6] = std::max(b[6], b[5]);
                    b[7] = std::max(b[7], b[5]);

                    b[4] = std::max(b[4], b[3]);
                    b[5] = std::max(b[5], b[3]);
                    b[6] = std::max(b[6], b[3]);
                    b[7] = std::max(b[7], b[3]);

                    a[i & 1][j + 0] = std::max(a[!(i & 1)][j + 0], b[0]);
                    a[i & 1][j + 1] = std::max(a[!(i & 1)][j + 1], b[1]);
                    a[i & 1][j + 2] = std::max(a[!(i & 1)][j + 2], b[2]);
                    a[i & 1][j + 3] = std::max(a[!(i & 1)][j + 3], b[3]);
                    a[i & 1][j + 4] = std::max(a[!(i & 1)][j + 4], b[4]);
                    a[i & 1][j + 5] = std::max(a[!(i & 1)][j + 5], b[5]);
                    a[i & 1][j + 6] = std::max(a[!(i & 1)][j + 6], b[6]);
                    a[i & 1][j + 7] = std::max(a[!(i & 1)][j + 7], b[7]);

                    j += 8;
                }

                for (; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), a[!(i & 1)][j - 1] + int(text1[i] == text2[j]));
                }
            }
        }
        return a[(text1.length() - 1) & 1][text2.length() - 1];
    }
};

// 4ms
class SolutionBy8_m128i : public Solution
{
private:
    int16_t a[5001][5001];
    int16_t b[8];

public:
    string name()
    {
        return "solution vectorized by 8";
    }

    int longestCommonSubsequence(string text1, string text2)
    {
        a[0][0] = int(text1[0] == text2[0]);
        // init entire first row and all columns
        for (int i = 1; i < text2.length(); i++)
            a[0][i] = std::max(a[0][i - 1], int16_t(text1[0] == text2[i]));

        if (text2.length() <= 8)
        {
            for (int i = 1; i < text1.length(); i++)
            {
                a[i & 1][0] = std::max(a[!(i & 1)][0], int16_t(text1[i] == text2[0]));
                for (int j = 1; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), int16_t(a[!(i & 1)][j - 1] + int16_t(text1[i] == text2[j])));
                }
            }
        }
        else
        {
            int16_t c[8];
            for (int i = 1; i < text1.length(); i++)
            {
                // get 8 values from previous row, c0 is dummy value?
                c[0] = int(text1[i] == text2[0]);
                c[1] = a[!(i & 1)][0] + int16_t(text1[i] == text2[1]);
                c[2] = a[!(i & 1)][1] + int16_t(text1[i] == text2[2]);
                c[3] = a[!(i & 1)][2] + int16_t(text1[i] == text2[3]);
                c[4] = a[!(i & 1)][3] + int16_t(text1[i] == text2[4]);
                c[5] = a[!(i & 1)][4] + int16_t(text1[i] == text2[5]);
                c[6] = a[!(i & 1)][5] + int16_t(text1[i] == text2[6]);
                c[7] = a[!(i & 1)][6] + int16_t(text1[i] == text2[7]);
                // load L(i - 1)(j - 1) + text1(i) == text2(i)
                __m128i b = _mm_loadu_si128((const __m128i_u *)c);

                // take max with L(i-1)(j)
                __m128i b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(13, 12, 13, 12, 9, 8, 9, 8, 5, 4, 5, 4, 1, 0, 1, 0));
                b = _mm_max_epi16(b, b_shifted);
                // b[1] = std::max(b[1], b[0]);
                // b[3] = std::max(b[3], b[2]);
                // b[5] = std::max(b[5], b[4]);
                // b[7] = std::max(b[7], b[6]);

                b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(11, 10, 11, 10, 11, 10, 9, 8, 3, 2, 3, 2, 3, 2, 1, 0));
                b = _mm_max_epi16(b, b_shifted);
                // b[2] = std::max(b[2], b[1]);
                // b[3] = std::max(b[3], b[1]);
                // b[6] = std::max(b[6], b[5]);
                // b[7] = std::max(b[7], b[5]);

                b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 5, 4, 3, 2, 1, 0));
                b = _mm_max_epi16(b, b_shifted);
                // b[4] = std::max(b[4], b[3]);
                // b[5] = std::max(b[5], b[3]);
                // b[6] = std::max(b[6], b[3]);
                // b[7] = std::max(b[7], b[3]);

                __m128i a0 = _mm_loadu_si128((const __m128i_u *)&a[!(i & 1)][0]);
                __m128i a1 = _mm_max_epi16(a0, b);
                _mm_storeu_si128((__m128i_u *)&a[i & 1][0], a1);
                // a[i&1][0] = std::max(a[!(i&1)][0], b[0]);
                // a[i&1][1] = std::max(a[!(i&1)][1], b[1]);
                // a[i&1][2] = std::max(a[!(i&1)][2], b[2]);
                // a[i&1][3] = std::max(a[!(i&1)][3], b[3]);
                // a[i&1][4] = std::max(a[!(i&1)][4], b[4]);
                // a[i&1][5] = std::max(a[!(i&1)][5], b[5]);
                // a[i&1][6] = std::max(a[!(i&1)][6], b[6]);
                // a[i&1][7] = std::max(a[!(i&1)][7], b[7]);

                int j = 8;
                while (j + 8 <= text2.length())
                {
                    int16_t k = _mm_extract_epi16(_mm_loadu_si128((__m128i *)&b), 7);
                    // int k = b[7];

                    __m128i s1 = _mm_set1_epi16(text1[i]);
                    __m128i s2 = _mm_loadu_si64(&text2[j]);
                    s2 = _mm_shuffle_epi8(s2, _mm_set_epi8(0x80, 7, 0x80, 6, 0x80, 5, 0x80, 4, 0x80, 3, 0x80, 2, 0x80, 1, 0x80, 0));
                    __m128i rs = _mm_cmpeq_epi16(s1, s2);
                    __m128i ones = _mm_set1_epi16(1);
                    rs = _mm_and_si128(rs, ones);
                    b = _mm_add_epi16(_mm_loadu_si128((__m128i *)&a[!(i & 1)][j - 1]), rs);
                    b = _mm_max_epi16(b, _mm_set1_epi16(k));
                    // b[0] = std::max(k, a[!(i&1)][j - 1] + int(text1[i]==text2[j]));
                    // b[1] = std::max(k, a[!(i&1)][j + 0] + int(text1[i]==text2[j + 1]));
                    // b[2] = std::max(k, a[!(i&1)][j + 1] + int(text1[i]==text2[j + 2]));
                    // b[3] = std::max(k, a[!(i&1)][j + 2] + int(text1[i]==text2[j + 3]));
                    // b[4] = std::max(k, a[!(i&1)][j + 3] + int(text1[i]==text2[j + 4]));
                    // b[5] = std::max(k, a[!(i&1)][j + 4] + int(text1[i]==text2[j + 5]));
                    // b[6] = std::max(k, a[!(i&1)][j + 5] + int(text1[i]==text2[j + 6]));
                    // b[7] = std::max(k, a[!(i&1)][j + 6] + int(text1[i]==text2[j + 7]));

                    __m128i b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(13, 12, 13, 12, 9, 8, 9, 8, 5, 4, 5, 4, 1, 0, 1, 0));
                    b = _mm_max_epi16(b, b_shifted);
                    // b[1] = std::max(b[1], b[0]);
                    // b[3] = std::max(b[3], b[2]);
                    // b[5] = std::max(b[5], b[4]);
                    // b[7] = std::max(b[7], b[6]);

                    b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(11, 10, 11, 10, 11, 10, 9, 8, 3, 2, 3, 2, 3, 2, 1, 0));
                    b = _mm_max_epi16(b, b_shifted);
                    // b[2] = std::max(b[2], b[1]);
                    // b[3] = std::max(b[3], b[1]);
                    // b[6] = std::max(b[6], b[5]);
                    // b[7] = std::max(b[7], b[5]);

                    b_shifted = _mm_shuffle_epi8(b, _mm_set_epi8(7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 5, 4, 3, 2, 1, 0));
                    b = _mm_max_epi16(b, b_shifted);
                    // b[4] = std::max(b[4], b[3]);
                    // b[5] = std::max(b[5], b[3]);
                    // b[6] = std::max(b[6], b[3]);
                    // b[7] = std::max(b[7], b[3]);

                    __m128i a0 = _mm_loadu_si128((const __m128i_u *)&a[!(i & 1)][j + 0]);
                    __m128i a1 = _mm_max_epi16(a0, b);
                    _mm_storeu_si128((__m128i_u *)&a[i & 1][j + 0], a1);
                    // a[i&1][j + 0] = std::max(a[!(i&1)][j + 0], b[0]);
                    // a[i&1][j + 1] = std::max(a[!(i&1)][j + 1], b[1]);
                    // a[i&1][j + 2] = std::max(a[!(i&1)][j + 2], b[2]);
                    // a[i&1][j + 3] = std::max(a[!(i&1)][j + 3], b[3]);
                    // a[i&1][j + 4] = std::max(a[!(i&1)][j + 4], b[4]);
                    // a[i&1][j + 5] = std::max(a[!(i&1)][j + 5], b[5]);
                    // a[i&1][j + 6] = std::max(a[!(i&1)][j + 6], b[6]);
                    // a[i&1][j + 7] = std::max(a[!(i&1)][j + 7], b[7]);

                    j += 8;
                }

                for (; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), int16_t(a[!(i & 1)][j - 1] + int16_t(text1[i] == text2[j])));
                }
            }
        }
        return a[(text1.length() - 1) & 1][text2.length() - 1];
    }
};

class SolutionBy8_m256i : public Solution
{
private:
    int16_t a[5001][5001];
    int16_t b[8];

public:
    string name()
    {
        return "solution vectorized by 8";
    }

    int longestCommonSubsequence(string text1, string text2)
    {
        a[0][0] = int(text1[0] == text2[0]);
        for (int i = 1; i < text2.length(); i++)
            a[0][i] = std::max(a[0][i - 1], int16_t(text1[0] == text2[i]));

        if (text2.length() <= 16)
        {
            for (int i = 1; i < text1.length(); i++)
            {
                a[i & 1][0] = std::max(a[!(i & 1)][0], int16_t(text1[i] == text2[0]));
                for (int j = 1; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), int16_t(a[!(i & 1)][j - 1] + int16_t(text1[i] == text2[j])));
                }
            }
        }
        else
        {
            for (int i = 1; i < text1.length(); i++)
            {
                // Load 16 characters from text2 and convert to 16-bit integers
                __m128i text2_bytes = _mm_loadu_si128((__m128i *)&text2[0]);
                __m256i text2_16bit = _mm256_cvtepu8_epi16(text2_bytes);
                // Broadcast text1[i] to all 16 positions
                __m256i text1_broadcast = _mm256_set1_epi16(text1[i]);
                // Compare for equality and convert to 0/1 values
                __m256i matches = _mm256_cmpeq_epi16(text1_broadcast, text2_16bit);
                __m256i ones = _mm256_set1_epi16(1);
                __m256i match_indicators = _mm256_and_si256(matches, ones);
                // Load two 128-bit chunks directly
                __m128i a_low = _mm_loadu_si128((__m128i *)&a[!(i & 1)][0]);  // a[0:7]
                __m128i a_high = _mm_loadu_si128((__m128i *)&a[!(i & 1)][8]); // a[8:15]
                // Shift right by one int16_t position: [a[0]...a[15]] -> [0, a[0]...a[14]]
                __m128i shifted_low = _mm_alignr_epi8(a_low, _mm_setzero_si128(), 14); // [0, a[0:6]]
                __m128i shifted_high = _mm_alignr_epi8(a_high, a_low, 14);             // [a[7:14]]
                __m256i a_shifted = _mm256_inserti128_si256(_mm256_castsi128_si256(shifted_low), shifted_high, 1);
                // Add match indicators to the shifted a values
                __m256i b = _mm256_add_epi16(match_indicators, a_shifted);
                // c[0] = int(text1[i]==text2[0]);
                // c[1] = a[!(i&1)][0] + int16_t(text1[i]==text2[1]);
                // c[2] = a[!(i&1)][1] + int16_t(text1[i]==text2[2]);
                // c[3] = a[!(i&1)][2] + int16_t(text1[i]==text2[3]);
                // c[4] = a[!(i&1)][3] + int16_t(text1[i]==text2[4]);
                // c[5] = a[!(i&1)][4] + int16_t(text1[i]==text2[5]);
                // c[6] = a[!(i&1)][5] + int16_t(text1[i]==text2[6]);
                // c[7] = a[!(i&1)][6] + int16_t(text1[i]==text2[7]);
                // c[8] = a[!(i&1)][7] + int16_t(text1[i]==text2[8]);
                // c[9] = a[!(i&1)][8] + int16_t(text1[i]==text2[9]);
                // c[10] = a[!(i&1)][9] + int16_t(text1[i]==text2[10]);
                // c[11] = a[!(i&1)][10] + int16_t(text1[i]==text2[11]);
                // c[12] = a[!(i&1)][11] + int16_t(text1[i]==text2[12]);
                // c[13] = a[!(i&1)][12] + int16_t(text1[i]==text2[13]);
                // c[14] = a[!(i&1)][13] + int16_t(text1[i]==text2[14]);
                // c[15] = a[!(i&1)][14] + int16_t(text1[i]==text2[15]);

                // First level: propagate every other element (0->1, 2->3, 4->5, etc.)
                __m256i b_shifted = _mm256_shuffle_epi8(b,
                                                        _mm256_set_epi8(29, 28, 29, 28, 25, 24, 25, 24, 21, 20, 21, 20, 17, 16, 17, 16,
                                                                        13, 12, 13, 12, 9, 8, 9, 8, 5, 4, 5, 4, 1, 0, 1, 0));
                b = _mm256_max_epi16(b, b_shifted);
                // b[1] = max(b[1], b[0]); b[3] = max(b[3], b[2]); etc.

                // Second level: propagate every 2 elements (0,1->2,3, 4,5->6,7, etc.)
                b_shifted = _mm256_shuffle_epi8(b,
                                                _mm256_set_epi8(27, 26, 27, 26, 27, 26, 25, 24, 19, 18, 19, 18, 19, 18, 17, 16,
                                                                11, 10, 11, 10, 11, 10, 9, 8, 3, 2, 3, 2, 3, 2, 1, 0));
                b = _mm256_max_epi16(b, b_shifted);
                // b[2] = max(b[2], b[1]); b[3] = max(b[3], b[1]); etc.

                // Third level: propagate every 4 elements (0,1,2,3->4,5,6,7, etc.)
                b_shifted = _mm256_shuffle_epi8(b,
                                                _mm256_set_epi8(23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 21, 20, 19, 18, 17, 16,
                                                                7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 5, 4, 3, 2, 1, 0));
                b = _mm256_max_epi16(b, b_shifted);
                // b[4] = max(b[4], b[3]); b[5] = max(b[5], b[3]); etc.

                // Fourth level: propagate every 8 elements (0-7->8-15)
                __m128i b_low = _mm256_castsi256_si128(b);
                __m128i b_high = _mm256_extracti128_si256(b, 1);
                __m128i max_low = _mm_set1_epi16(_mm_extract_epi16(b_low, 7));
                b_high = _mm_max_epi16(b_high, max_low);
                b = _mm256_inserti128_si256(b, b_high, 1);

                __m256i a0 = _mm256_loadu_si256((const __m256i *)&a[!(i & 1)][0]);
                __m256i a1 = _mm256_max_epi16(a0, b);
                _mm256_storeu_si256((__m256i *)&a[i & 1][0], a1);

                int j = 16;
                while (j + 16 <= text2.length())
                {
                    int16_t k = _mm_extract_epi16(_mm256_extracti128_si256(b, 1), 7);
                    // int k = b[15]; // highest element from previous iteration

                    // Load 16 characters from text2 and compare with text1[i]
                    __m256i s1 = _mm256_set1_epi16(text1[i]);
                    __m128i s2_bytes = _mm_loadu_si128((__m128i *)&text2[j]);
                    __m256i s2 = _mm256_cvtepu8_epi16(s2_bytes);
                    __m256i rs = _mm256_cmpeq_epi16(s1, s2);
                    __m256i ones = _mm256_set1_epi16(1);
                    rs = _mm256_and_si256(rs, ones);
                    b = _mm256_add_epi16(_mm256_loadu_si256((__m256i *)&a[!(i & 1)][j - 1]), rs);
                    b = _mm256_max_epi16(b, _mm256_set1_epi16(k));

                    // Apply the same propagation pattern as before
                    __m256i b_shifted = _mm256_shuffle_epi8(b,
                                                            _mm256_set_epi8(29, 28, 29, 28, 25, 24, 25, 24, 21, 20, 21, 20, 17, 16, 17, 16,
                                                                            13, 12, 13, 12, 9, 8, 9, 8, 5, 4, 5, 4, 1, 0, 1, 0));
                    b = _mm256_max_epi16(b, b_shifted);

                    b_shifted = _mm256_shuffle_epi8(b,
                                                    _mm256_set_epi8(27, 26, 27, 26, 27, 26, 25, 24, 19, 18, 19, 18, 19, 18, 17, 16,
                                                                    11, 10, 11, 10, 11, 10, 9, 8, 3, 2, 3, 2, 3, 2, 1, 0));
                    b = _mm256_max_epi16(b, b_shifted);

                    b_shifted = _mm256_shuffle_epi8(b,
                                                    _mm256_set_epi8(23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 21, 20, 19, 18, 17, 16,
                                                                    7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 5, 4, 3, 2, 1, 0));
                    b = _mm256_max_epi16(b, b_shifted);

                    // Fourth level propagation for 256-bit
                    __m128i b_low = _mm256_castsi256_si128(b);
                    __m128i b_high = _mm256_extracti128_si256(b, 1);
                    __m128i max_low = _mm_set1_epi16(_mm_extract_epi16(b_low, 7));
                    b_high = _mm_max_epi16(b_high, max_low);
                    b = _mm256_inserti128_si256(b, b_high, 1);

                    __m256i a0 = _mm256_loadu_si256((const __m256i *)&a[!(i & 1)][j]);
                    __m256i a1 = _mm256_max_epi16(a0, b);
                    _mm256_storeu_si256((__m256i *)&a[i & 1][j], a1);

                    j += 16;
                }

                for (; j < text2.length(); j++)
                {
                    a[i & 1][j] = std::max(std::max(a[!(i & 1)][j], a[i & 1][j - 1]), int16_t(a[!(i & 1)][j - 1] + int16_t(text1[i] == text2[j])));
                }
            }
        }
        return a[(text1.length() - 1) & 1][text2.length() - 1];
    }
};

void benchmark(Solution *sol, int n, int repeat_time = 100)
{
    string text1;
    string text2;
    for (int i = 0; i < n; i++)
    {
        text1.push_back('a' + i % 26);
        text2.push_back('a' + i % 25);
    }

    sol->longestCommonSubsequence(text1, text2);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < repeat_time; i++)
    {
        sol->longestCommonSubsequence(text1, text2);
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << sol->name() << " execution time: " << duration.count() << " seconds\n";
}

int main()
{
    Solution *sol = new SolutionNaive();
    benchmark(sol, 1024);
    sol = new SolutionCompress();
    benchmark(sol, 1024);
    sol = new SolutionBy4();
    benchmark(sol, 1024);
    sol = new SolutionBy8();
    benchmark(sol, 1024);
    sol = new SolutionBy8_m256i();
    benchmark(sol, 1024);
    sol = new SolutionBy8_m128i();
    benchmark(sol, 1024);
    return 0;
}