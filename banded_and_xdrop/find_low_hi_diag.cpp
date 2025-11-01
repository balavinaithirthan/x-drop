#include <vector>
#include <string>
#include <algorithm>
#include <limits>

using Matrix = std::vector<std::vector<double>>;
constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double POS_INF =  std::numeric_limits<double>::infinity();

Matrix nw4(const std::string& rowSeq,
           const std::string& colSeq,
           int rows, int cols,
           Matrix mat,
           double x_thresh)
{
    rows += 1;
    cols += 1;

    int i_start = 0;
    int i_end   = 0;
    double max_score  = NEG_INF;
    double lower_diag = NEG_INF;
    double upper_diag = POS_INF;

    for (int ad = 0; ad < rows + cols - 1; ++ad) {
        i_end   = std::min(i_end + 1, rows);
        i_start = std::max(i_start, ad - cols + 1);
        if (i_start >= i_end) break;

        for (int i = i_start; i < i_end; ++i) {
            int j = ad - i;

            double diag = mat[i - 1][j - 1] + (rowSeq[i - 1] == colSeq[j - 1] ? 1 : 0);
            double left = mat[i][j - 1] - 2;
            double up   = mat[i - 1][j] - 2;
            double s    = std::max({diag, left, up});
            mat[i][j] = s;
            max_score = std::max(max_score, s);
        }

        for (int i = i_start; i < i_end; ++i) {
            int j = ad - i;
            if (mat[i][j] < max_score - x_thresh)
                mat[i][j] = NEG_INF;
        }

        // get_i_start logic
        {
            int i = i_start;
            int j = ad - i;
            while (i < i_end && mat[i][j] == NEG_INF) {
                ++i;
                j = ad - i;
            }
            double d_lo;
            if (mat[i_start][ad - i_start] != NEG_INF)
                d_lo = NEG_INF;
            else
                d_lo = i - j;
            i_start = i;
            lower_diag = std::max(lower_diag, d_lo);
        }

        // get_i_end logic
        {
            int i = i_end - 1;
            int j = ad - i;
            while (i > i_start && mat[i][j] == NEG_INF) {
                --i;
                j = ad - i;
            }
            double d_hi;
            if (mat[i_end - 1][ad - (i_end - 1)] != NEG_INF)
                d_hi = POS_INF;
            else
                d_hi = i - j + 2;
            i_end = i + 1;
            upper_diag = std::min(upper_diag, d_hi);
        }

        // modify_i_start logic
        {
            int i_on_diag = static_cast<int>((lower_diag + ad) / 2);
            i_start = std::max(i_start, i_on_diag);
        }

        // modify_i_end logic
        {
            int i_on_diag = static_cast<int>((upper_diag + ad) / 2);
            i_end = std::min(i_end, i_on_diag);
        }
    }

    return mat;
}
