#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <iomanip> // for std::setw
#include <random>

// We'll use a global 'ping' like the Python code.
static const int ping = 1;

/**
 * Initialize an N x N matrix (where N = len(A) + ping) with 'fill_val' in its
 * interior, but sets the base-case boundaries: mat[i][0] = -i, mat[0][j] = -j.
 *
 * NOTE: In the original Python code, for demonstration, the matrix is sized
 * by A’s length alone, ignoring B’s length. This works for the example if A
 * and B are the same length. Adjust for the general case if needed.
 */
std::vector<std::vector<double>> initializeMatrix(const std::string &A,
                                                  const std::string &B,
                                                  double fill_val)
{
    int N = static_cast<int>(A.size());
    // We'll replicate the original Python approach exactly:
    std::vector<std::vector<double>> mat(N + ping, std::vector<double>(N + ping, 0.0));

    // Fill interior (excluding the top row & left col) with fill_val
    // the Python code does matp = mat[1:,1:], matp.fill(fill_val)
    for (int i = 1; i < N + ping; i++)
    {
        for (int j = 1; j < N + ping; j++)
        {
            mat[i][j] = fill_val;
        }
    }

    // Base case
    for (int i = 1; i < static_cast<int>(A.size()) + ping; i++)
    {
        mat[i][0] = -static_cast<double>(i);
    }
    for (int j = 1; j < static_cast<int>(B.size()) + ping; j++)
    {
        // Careful not to exceed the matrix boundary if A.size() != B.size().
        if (j < N + ping)
            mat[0][j] = -static_cast<double>(j);
    }
    return mat;
}

/**
 * Print a matrix for debugging.
 */
void printMatrix(const std::vector<std::vector<double>> &mat)
{
    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < mat[i].size(); j++)
        {
            std::cout << std::setw(6) << mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/**
 * Score function:
 *   diag = mat[i-1][j-1] + (1 if A[i-1] == B[j-1], else 0)
 *   left = mat[i][j-1] - 2
 *   up   = mat[i-1][j] - 2
 * Return max(diag, left, up).
 *
 * We assume i and j are within valid index range in the matrix.
 */
double score(const std::vector<std::vector<double>> &mat,
             const std::string &A,
             const std::string &B,
             int i,
             int j)
{
    double diag = mat[i - 1][j - 1] + (A[i - 1] == B[j - 1] ? 1.0 : 0.0);
    double left = mat[i][j - 1] - 2.0;
    double up = mat[i - 1][j] - 2.0;
    return std::max({diag, left, up});
}

/**
 * get_i_start() from the Python code.
 *   We move i upward from i_start while mat[i][j] == -∞, then return the new i
 *   and a diagonal offset d_lo.
 */
std::pair<int, double> get_i_start(
    std::vector<std::vector<double>> &mat,
    std::vector<std::vector<double>> &explored,
    int d,
    int i_start,
    int i_end)
{
    int i = i_start;
    double d_lo = -std::numeric_limits<double>::infinity();
    while (i < i_end)
    {
        int j = d - i;
        if (j < 0 || j >= static_cast<int>(mat[i].size()))
            break;
        if (mat[i][j] != -std::numeric_limits<double>::infinity())
        {
            // We found a valid cell
            break;
        }
        // Mark explored
        explored[i][j] = 1.0;
        i++;
    }
    // If the original i_start cell is not -∞, set d_lo to -∞
    // else set something else. (Mirroring the Python code logic.)
    if (i_start < i_end)
    {
        int j_start = d - i_start;
        if (0 <= j_start && j_start < (int)mat[i_start].size())
        {
            if (mat[i_start][j_start] != -std::numeric_limits<double>::infinity())
            {
                d_lo = -std::numeric_limits<double>::infinity();
            }
            else
            {
                // The Python code does: d_lo = i - j + ...
                // "i - j" is "i - (d - i)" = i - d + i = 2*i - d
                // They used +2 in some place. We'll replicate carefully:
                int j_temp = d - i;
                d_lo = (double)(2 * i - d);
            }
        }
    }
    return std::make_pair(i, d_lo);
}

/**
 * get_i_end() from Python code. We move downward from i_end-1
 * while mat[i][j] == -∞, then return new i_end and d_hi.
 */
std::pair<int, double> get_i_end(
    std::vector<std::vector<double>> &mat,
    std::vector<std::vector<double>> &explored,
    int d,
    int i_start,
    int i_end)
{
    int i = i_end - 1; // we walk downward from i_end-1
    double d_hi = std::numeric_limits<double>::infinity();
    while (i >= i_start)
    {
        int j = d - i;
        if (j < 0 || j >= static_cast<int>(mat[i].size()))
            break;
        if (mat[i][j] != -std::numeric_limits<double>::infinity())
        {
            // Found a valid cell
            break;
        }
        explored[i][j] = 1.0;
        i--;
    }
    // If mat[i_end-1][?] != -∞, set d_hi = ∞ else something else
    int i_chk = i_end - 1;
    int j_chk = d - i_chk;
    if (i_chk >= 0 && i_chk < (int)mat.size() && j_chk >= 0 && j_chk < (int)mat[i_chk].size())
    {
        if (mat[i_chk][j_chk] != -std::numeric_limits<double>::infinity())
        {
            d_hi = std::numeric_limits<double>::infinity();
        }
        else
        {
            // The python code: d_hi = i - j + 2 => i - (d - i) + 2 = 2*i - d + 2
            int j_temp = d - i;
            d_hi = (double)(2 * i - d + 2);
        }
    }
    return std::make_pair(i + 1, d_hi);
}

/**
 * Mark cells whose score < (max_score - X_drop) as -∞ and also mark 'explored'.
 */
void markX(std::vector<std::vector<double>> &mat,
           std::vector<std::vector<double>> &explored,
           int d, int i_start, int i_end,
           double X_drop,
           double max_score)
{
    // replicate the "for i in range(i_start, i_end):" logic
    for (int i = i_start; i < i_end; i++)
    {
        int j = d - i;
        if (j < 0 || j >= (int)mat[i].size())
            continue;
        double s = mat[i][j];
        explored[i][j] = 1.0;
        if (s < (max_score - X_drop))
        {
            mat[i][j] = -std::numeric_limits<double>::infinity();
        }
    }
}

/**
 * modify_i_start: from python code
 *   i_on_diag = (d_lo + ad) // 2
 *   i_start = max(i_start, i_on_diag)
 */
int modify_i_start(int i_start, double d_lo, int ad)
{
    if (d_lo == -std::numeric_limits<double>::infinity())
    {
        // same as python: if no real offset, just keep i_start
        return i_start;
    }
    // The python code does (d_lo + ad) // 2 for i_on_diag
    // We'll do integer division in C++ as well:
    int i_on_diag = (int)((d_lo + (double)ad) / 2.0);
    return std::max(i_start, i_on_diag);
}

/**
 * modify_i_end: from python code
 *   i_on_diag = (d_hi + ad) // 2
 *   i_end = min(i_end, i_on_diag)
 */
int modify_i_end(int i_end, double d_hi, int ad)
{
    if (d_hi == std::numeric_limits<double>::infinity())
    {
        return i_end;
    }
    int i_on_diag = (int)((d_hi + (double)ad) / 2.0);
    return std::min(i_end, i_on_diag);
}

/**
 * The main x-drop NW function, replicating the Python's nw4(A, B, x_thresh=2).
 * Returns the final matrix with computed scores and the explored matrix.
 *
 * The logic is specialized for the example, but follows the Python code
 * anti-diagonal iteration.
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> xdrop_nw(const std::string &A,
                                                                                       const std::string &B,
                                                                                       double x_thresh)
{

    // Matrix dimension (following the python snippet):
    std::vector<std::vector<double>> mat = initializeMatrix(A, B, 0.1);

    // We'll also keep an 'explored' matrix to track visited cells
    std::vector<std::vector<double>> explored = initializeMatrix(A, B, 0.0);

    int m = static_cast<int>(A.size()) + 1; // rows
    int n = static_cast<int>(B.size()) + 1; // cols

    int i_start = 1;
    int i_end = 1;
    double max_score = -std::numeric_limits<double>::infinity();

    double lower_diag = -std::numeric_limits<double>::infinity();
    double upper_diag = std::numeric_limits<double>::infinity();
    int final_ad = 0;
    // First pass: ad in [2, m)
    for (int ad = 2; ad < m; ad++)
    {
        i_end++;
        for (int i = i_start; i < i_end; i++)
        {
            int j = ad - i;
            if (j < 0 || j >= (int)mat[i].size())
            {
                continue;
            }
            // compute the score
            double s = score(mat, A, B, i, j);
            mat[i][j] = s;
            explored[i][j] = 1.0;
            if (s > max_score)
            {
                max_score = s;
            }
        }
        // Mark X-drop
        markX(mat, explored, ad, i_start, i_end, x_thresh, max_score);

        // Adjust i_start, i_end
        auto start_pair = get_i_start(mat, explored, ad, i_start, i_end);
        i_start = start_pair.first;
        double d_lo = start_pair.second;

        auto end_pair = get_i_end(mat, explored, ad, i_start, i_end);
        i_end = end_pair.first;
        double d_hi = end_pair.second;

        lower_diag = std::max(lower_diag, d_lo);
        upper_diag = std::min(upper_diag, d_hi);

        i_start = modify_i_start(i_start, lower_diag, ad);
        i_end = modify_i_end(i_end, upper_diag, ad);
    }

    // Second pass: ad in [m, m + n - 1)
    // This covers anti-diagonals that extend beyond the first dimension's corner.
    for (int ad = m; ad < m + n - 1; ad++)
    {
        int expected_i_start = ad - m + 1;
        if (expected_i_start > i_start)
        {
            i_start = expected_i_start;
        }
        int expected_i_end = m;
        if ((i_end + 1) < expected_i_end)
        {
            i_end = i_end + 1;
        }
        else
        {
            i_end = expected_i_end;
        }

        for (int i = i_start; i < i_end; i++)
        {
            int j = ad - i;
            if (j < 0 || j >= (int)mat[i].size() || i < 0 || i >= (int)mat.size())
            {
                continue;
            }
            double s = score(mat, A, B, i, j);
            mat[i][j] = s;
            explored[i][j] = 1.0;
            if (s > max_score)
            {
                max_score = s;
            }
        }

        if (i_start == i_end)
        {
            // no more space to move -> break
            break;
        }

        // Mark X-drop
        markX(mat, explored, ad, i_start, i_end, x_thresh, max_score);

        // Adjust i_start, i_end
        auto start_pair = get_i_start(mat, explored, ad, i_start, i_end);
        i_start = start_pair.first;
        double d_lo = start_pair.second;

        auto end_pair = get_i_end(mat, explored, ad, i_start, i_end);
        i_end = end_pair.first;
        double d_hi = end_pair.second;

        lower_diag = std::max(lower_diag, d_lo);
        upper_diag = std::min(upper_diag, d_hi);

        i_start = modify_i_start(i_start, lower_diag, ad);
        i_end = modify_i_end(i_end, upper_diag, ad);
        final_ad = ad;
    }
    printf("final antidiagonal explored is %d\n", final_ad);
    return std::make_pair(mat, explored);
}

std::string string_gen(const int length)
{
    // The alphabet of possible characters
    static const char alphabet[] = {'A', 'G', 'C', 'T'};

    // Set up random number generation
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister
    std::uniform_int_distribution<> dist(0, 3);

    // Build a random string
    std::string randomDNA;
    randomDNA.reserve(length); // reserve for efficiency

    for (int i = 0; i < length; ++i)
    {
        randomDNA.push_back(alphabet[dist(gen)]);
    }
    return randomDNA;
}

int main()
{
    std::string A = string_gen(10000);
    std::string B = string_gen(10000);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto result = xdrop_nw(A, B, 2.0);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    std::vector<std::vector<double>> mat = result.first;
    std::vector<std::vector<double>> explored = result.second;

    // Print out the final matrix
    std::cout << "Final x-drop NW matrix:\n";
    std::cout << mat[A.size()][B.size()] << std::endl;

    // Count the number of explored cells where the value is 1.0
    int explored_count = 0;
    for (const auto &row : explored)
    {
        for (const auto &cell : row)
        {
            if (cell == 1.0)
            {
                explored_count++;
            }
        }
    }

    std::cout << "Number of explored cells: " << explored_count << std::endl;
    std::cout << "Totoal cells in the matrix: " << A.size() * B.size() << std::endl;
    return 0;
}
