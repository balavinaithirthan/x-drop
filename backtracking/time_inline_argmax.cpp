#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

struct AlignmentResult
{
    int max_score;
    int max_i;
    int max_j;
    std::vector<std::vector<int>> matrix;
};

// Smith-Waterman with inline maximum tracking
AlignmentResult smith_waterman_inline_max(std::vector<std::vector<int>> &H, const std::string &seq1, const std::string &seq2,
                                          int match = 2, int mismatch = -1, int gap = -1)
{
    int m = seq1.length();
    int n = seq2.length();

    // Initialize scoring matrix

    // Track maximum score and position inline
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;
    int count_cells = 0;

    // Fill the scoring matrix
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // Calculate scores for each possible move
            int diagonal = H[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int up = H[i - 1][j] + gap;
            int left = H[i][j - 1] + gap;

            // Take maximum of all possibilities, including 0 (local alignment)
            H[i][j] = std::max({0, diagonal, up, left});
            count_cells++;

            // Track maximum score and position inline (hot loop)
            if (H[i][j] > max_score)
            {
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    return {max_score, max_i, max_j, H};
}

// Smith-Waterman with post-computation maximum finding
AlignmentResult smith_waterman_post_max(std::vector<std::vector<int>> &H, const std::string &seq1, const std::string &seq2,
                                        int match = 2, int mismatch = -1, int gap = -1)
{
    int m = seq1.length();
    int n = seq2.length();
    int count_cells = 0;
    // Initialize scoring matrix

    // Fill the scoring matrix (hot loop without maximum tracking)
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // Calculate scores for each possible move
            int diagonal = H[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int up = H[i - 1][j] + gap;
            int left = H[i][j - 1] + gap;

            // Take maximum of all possibilities, including 0 (local alignment)
            H[i][j] = std::max({0, diagonal, up, left});
            count_cells++;
        }
    }

    // Find maximum score and position after matrix computation
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            if (H[i][j] > max_score)
            {
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    return {max_score, max_i, max_j, H};
}

// Function to print the alignment matrix
void print_matrix(const std::vector<std::vector<int>> &matrix, const std::string &seq1, const std::string &seq2)
{
    std::cout << "    ";
    for (char c : seq2)
    {
        std::cout << c << "  ";
    }
    std::cout << std::endl;

    for (int i = 0; i < matrix.size(); i++)
    {
        if (i == 0)
        {
            std::cout << "  ";
        }
        else
        {
            std::cout << seq1[i - 1] << " ";
        }

        for (int j = 0; j < matrix[i].size(); j++)
        {
            std::cout << matrix[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}

// Function to benchmark the two approaches
void benchmark_functions(const std::string &seq1, const std::string &seq2, int iterations = 1000)
{
    std::cout << "Benchmarking Smith-Waterman implementations..." << std::endl;
    std::cout << "Sequence 1 length: " << seq1.length() << std::endl;
    std::cout << "Sequence 2 length: " << seq2.length() << std::endl;
    std::cout << "Iterations: " << iterations << std::endl
              << std::endl;

    // Benchmark inline max tracking
    std::vector<std::vector<int>> H = std::vector<std::vector<int>>(seq1.length() + 1, std::vector<int>(seq2.length() + 1, 0));
    auto start = std::chrono::high_resolution_clock::now();
    AlignmentResult result_inline;
    for (int i = 0; i < iterations; i++)
    {
        result_inline = smith_waterman_inline_max(H, seq1, seq2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_inline = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Benchmark post-computation max finding
    H = std::vector<std::vector<int>>(seq1.length() + 1, std::vector<int>(seq2.length() + 1, 0));

    start = std::chrono::high_resolution_clock::now();
    AlignmentResult result_post;
    for (int i = 0; i < iterations; i++)
    {
        result_post = smith_waterman_post_max(H, seq1, seq2);
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration_post = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Print results
    std::cout << "Results:" << std::endl;
    std::cout << "Inline max tracking - Max score: " << result_inline.max_score
              << " at position (" << result_inline.max_i << ", " << result_inline.max_j << ")" << std::endl;
    std::cout << "Post-computation max - Max score: " << result_post.max_score
              << " at position (" << result_post.max_i << ", " << result_post.max_j << ")" << std::endl;

    std::cout << std::endl
              << "Performance:" << std::endl;
    std::cout << "Inline max tracking: " << duration_inline.count() << " microseconds" << std::endl;
    std::cout << "Post-computation max: " << duration_post.count() << " microseconds" << std::endl;
    std::cout << "Speedup: " << (double)duration_post.count() / duration_inline.count() << "x" << std::endl;
}

int main()
{
    // Test sequences
    std::string seq1 = std::string("A", 3000); // 10000 long sequence
    std::string seq2 = std::string("B", 3000); // 10000 long sequence

    // std::cout << "Testing Smith-Waterman implementations" << std::endl;
    // std::cout << "=======================================" << std::endl
    //           << std::endl;

    // // Test with identical sequences
    // std::cout << "Test 1: Identical sequences" << std::endl;
    // std::cout << "Sequence 1: " << seq1 << std::endl;
    // std::cout << "Sequence 2: " << seq2 << std::endl
    //           << std::endl;

    // AlignmentResult result1 = smith_waterman_inline_max(seq1, seq2);
    // AlignmentResult result2 = smith_waterman_post_max(seq1, seq2);

    // std::cout << "Inline max - Score: " << result1.max_score
    //           << " at (" << result1.max_i << ", " << result1.max_j << ")" << std::endl;
    // std::cout << "Post max - Score: " << result2.max_score
    //           << " at (" << result2.max_i << ", " << result2.max_j << ")" << std::endl;

    // std::cout << std::endl
    //           << "Matrix (inline max):" << std::endl;
    // print_matrix(result1.matrix, seq1, seq2);

    // // Test with different sequences
    // std::cout << std::endl
    //           << "Test 2: Different sequences" << std::endl;
    // seq1 = "GGTTGACTA";
    // seq2 = "TGTTACGG";
    // std::cout << "Sequence 1: " << seq1 << std::endl;
    // std::cout << "Sequence 2: " << seq2 << std::endl
    //           << std::endl;

    // result1 = smith_waterman_inline_max(seq1, seq2);
    // result2 = smith_waterman_post_max(seq1, seq2);

    // std::cout << "Inline max - Score: " << result1.max_score
    //           << " at (" << result1.max_i << ", " << result1.max_j << ")" << std::endl;
    // std::cout << "Post max - Score: " << result2.max_score
    //           << " at (" << result2.max_i << ", " << result2.max_j << ")" << std::endl;

    // std::cout << std::endl
    //           << "Matrix (inline max):" << std::endl;
    // print_matrix(result1.matrix, seq1, seq2);

    // Benchmark performance
    std::cout << std::endl
              << "Performance Benchmark" << std::endl;
    std::cout << "====================" << std::endl;
    benchmark_functions(seq1, seq2, 10);

    return 0;
}
