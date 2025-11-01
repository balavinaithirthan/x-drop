#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <random>
#include <chrono>

// Score parameters for NW-like scoring
static constexpr double MATCH_SCORE   =  1.0;
static constexpr double MISMATCH_PEN  =  0.0;  // (You could set e.g. -1.0 if needed)
static constexpr double GAP_PENALTY   = -2.0;

// For convenience, represent -∞ as:
static const double NEG_INF = -std::numeric_limits<double>::infinity();

/**
 * Return the alignment score of A[i-1], B[j-1].
 * In a simple match/mismatch scheme:
 *     match => +1
 *     mismatch => +0 (or negative, depending on your scoring)
 */
inline double matchScore(char a, char b)
{
    return (a == b) ? MATCH_SCORE : MISMATCH_PEN;
}

/**
 * Generate a random DNA string of given length.
 */
std::string string_gen(int length)
{
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 3);

    std::string s; 
    s.reserve(length);
    for(int i = 0; i < length; i++)
        s.push_back(alphabet[dist(gen)]);
    return s;
}

/**
 * x-drop NW using only three anti-diagonals in memory.
 *
 * We treat d = i + j as the “index” of the current anti-diagonal. For each d:
 *   - i ranges between i_start and i_end
 *   - j = d - i
 * We store:
 *   diag(d)   in a 1D vector named curr_diag
 *   diag(d-1) in a 1D vector named prev_diag
 *   diag(d-2) in a 1D vector named prev2_diag
 *
 * The code performs an x-drop: any cell whose score < (max_score - x_drop)
 * is pruned (set to -∞), and future expansions from that cell are cut off.
 *
 * Return value: the best (max) alignment score found, as well as
 * the total number of DP cells “visited” (for example, to measure complexity).
 */
std::pair<double, long long> xdrop_nw_3diag(const std::string &A,
                                            const std::string &B,
                                            double x_drop)
{
    // We'll define these for convenience:
    int lenA = (int)A.size();
    int lenB = (int)B.size();

    // We want to handle boundaries: DP[0][j] = -j, DP[i][0] = -i.
    // We'll define “d” to go from 0 up to (lenA + lenB).
    // For each diagonal d, i + j = d => j = d - i.
    // i goes from max(0, d-lenB) to min(d, lenA).

    // Keep track of the best score so far:
    double max_score = NEG_INF;

    // We will store three diagonals:
    //   prev2_diag => diag(d-2)
    //   prev_diag  => diag(d-1)
    //   curr_diag  => diag(d)
    // The maximum possible diagonal length is min(lenA+1, lenB+1),
    // but for x-drop we’ll typically keep a narrower band.

    std::vector<double> prev2_diag, prev_diag, curr_diag;

    // i_start, i_end define the range of i for the current diagonal
    // We start at d=0 => DP[0][0].
    // Over time we expand i_end, then eventually shrink it as we pass
    // the corners. We'll track these to mimic your original code's approach.
    int i_start = 0;
    int i_end   = 1;  // for d=0, we only have (i=0, j=0)

    // We'll also keep track of "lower_diag" and "upper_diag" offsets
    // for x-drop bounding. (These replicate your original logic that
    // tries to clamp i_start and i_end.)
    double lower_diag = NEG_INF;
    double upper_diag =  std::numeric_limits<double>::infinity();

    // We'll keep track of how many cells we actually compute:
    long long explored_count = 0;

    // A small helper to resize a diagonal array for the new length
    // and fill with -∞.
    auto resize_fill_neg_inf = [&](std::vector<double> &diag, int length)
    {
        diag.assign(length, NEG_INF);
    };

    // A helper to compute DP[i][j] from the three neighbors:
    auto compute_cell_score = [&](int i, int d,
                                  const std::vector<double> &prev2,
                                  const std::vector<double> &prev,
                                  std::vector<double> &curr,
                                  int i_start_prev2,
                                  int i_start_prev,
                                  int i_start_curr)
    {
        // j = d - i
        int j = d - i;

        // The index of DP(i,j) in the current diagonal array:
        int idx_curr = i - i_start_curr;

        // diag (i-1, j-1) => from diagonal d-2
        // i-1 => index in prev2 = (i-1) - i_start_prev2
        double diag_val = NEG_INF;
        if(i > 0 && j > 0) {
            int idx_diag = (i - 1) - i_start_prev2;
            if(idx_diag >= 0 && idx_diag < (int)prev2.size()) {
                diag_val = prev2[idx_diag];
                if(diag_val != NEG_INF) {
                    diag_val += matchScore(A[i-1], B[j-1]);
                }
            }
        }

        // left => from diagonal d-1 at the same i => DP[i][j-1]
        double left_val = NEG_INF;
        if(j > 0) {
            int idx_left = i - i_start_prev;
            if(idx_left >= 0 && idx_left < (int)prev.size()) {
                left_val = prev[idx_left];
                if(left_val != NEG_INF) {
                    left_val += GAP_PENALTY; // j-1 => gap in B
                }
            }
        }

        // up => from diagonal d-1 at i-1 => DP[i-1][j]
        double up_val = NEG_INF;
        if(i > 0) {
            int idx_up = (i - 1) - i_start_prev;
            if(idx_up >= 0 && idx_up < (int)prev.size()) {
                up_val = prev[idx_up];
                if(up_val != NEG_INF) {
                    up_val += GAP_PENALTY; // i-1 => gap in A
                }
            }
        }

        double best = std::max({diag_val, left_val, up_val});
        curr[idx_curr] = best;
        return best;
    };

    // We also need to replicate the boundary condition:
    //   DP[i][0] = -i,  DP[0][j] = -j.
    // For an anti‐diagonal approach, that means:
    //  - When d=0, DP[0][0] = 0
    //  - On future diagonals, if i=0 => DP[0][j] = -j
    //    if j=0 => DP[i][0] = -i
    auto boundary_condition = [&](int i, int j){
        if(i == 0 && j >= 0) return -double(j);
        if(j == 0 && i >= 0) return -double(i);
        return NEG_INF; 
    };

    // For convenience, define small lambdas that replicate the
    // "get_i_start" / "get_i_end" and "mark X-drop" logic from your code.
    // We'll do simpler versions: 
    //   - We'll skip the complicated 'd_lo' / 'd_hi' usage,
    //   - We'll do a linear sweep from i_start to i_end to see if the front or back
    //     are -∞ for bounding.
    // (You can adapt the more complicated logic as needed.)

    auto trim_i_start = [&](std::vector<double> &diag, int &start, int &end) {
        // Move 'start' up while diag is -∞
        while(start < end) {
            int idx = start - i_start; // but we must know how we set i_start for diag
            if(idx < 0 || idx >= (int)diag.size()) break;
            if(diag[idx] != NEG_INF) break; // found a valid cell
            start++;
        }
    };

    auto trim_i_end = [&](std::vector<double> &diag, int &start, int &end) {
        // Move 'end' down while diag is -∞
        while(end > start) {
            int idx = end - 1 - i_start; 
            if(idx < 0 || idx >= (int)diag.size()) break;
            if(diag[idx] != NEG_INF) break; // found a valid cell
            end--;
        }
    };

    // For x-drop pruning: we mark any cell < (max_score - x_drop) as -∞
    auto xdrop_prune = [&](std::vector<double> &diag) {
        for(auto &val : diag) {
            if(val < max_score - x_drop) {
                val = NEG_INF;
            }
        }
    };

    // ---- Main iteration over diagonals ----
    double best_score = NEG_INF;

    // We will handle each diagonal d from 0 to (lenA + lenB).
    // We do it in two phases as your original code does, or we can do it in one pass.
    // For simplicity, we'll just do it in one pass from d = 0..(lenA + lenB).

    // Initialize first diagonal, d=0 => DP[0][0] = 0
    prev2_diag.clear();
    prev_diag.clear();
    curr_diag.clear();
    resize_fill_neg_inf(curr_diag, 1);
    curr_diag[0] = 0.0;  // DP[0][0] = 0
    max_score = 0.0;
    explored_count = 1; // we computed one cell
    best_score = 0.0;

    // We will define a few helper variables that track i_start for each diagonal
    int i_start_prev2 = 0, i_start_prev = 0;
    i_start = 0; 
    // for d=0, i_end=1

    // now shift diagonals: prev2 <- empty, prev <- curr
    prev2_diag = std::move(prev2_diag); // empty
    prev_diag  = std::move(curr_diag);

    i_start_prev2 = -1; // meaningless for now
    i_start_prev  = 0;
    int final_ad = 0;
    // We proceed from d=1 up to d = lenA + lenB
    for(int d = 1; d <= (lenA + lenB); d++)
    {
        // The range of i is [max(0, d - lenB), min(d, lenA)]
        // but we also incorporate any running i_start/i_end from x-drop bounding
        int new_i_start = std::max(0, d - lenB);
        int new_i_end   = std::min(d, lenA);

        // Possibly expand or shrink based on the prior iteration
        // In a standard NW, you keep the full diagonal.  For x-drop,
        // you might skip portions.  For simplicity, we'll just do:
        //   i_start = std::max(i_start, new_i_start);
        //   i_end   = std::min(i_end+1, new_i_end);
        // or replicate the two-phase logic. This is a simplified approach:

        if(d <= lenA) {
            // expanding the diagonal downward
            if(d == 1) {
                i_start = 0;
                i_end   = 2; // from the previous pass, i_end was 1
            } else {
                i_end++;
            }
            if(i_end > new_i_end + 1) {
                i_end = new_i_end + 1;
            }
        } else {
            // after we've passed the 'corner', we fix i_end at lenA+1
            // and eventually shrink from the top if d-lenB grows
            if(i_end < new_i_end + 1) {
                i_end++;
            }
            if(i_end > (new_i_end + 1)) {
                i_end = new_i_end + 1;
            }
        }

        if(i_start < new_i_start) {
            i_start = new_i_start;
        }
        if(i_end > new_i_end + 1) {
            i_end = new_i_end + 1;
        }
        if(i_end < i_start) {
            // no valid cells, we can break early
            break;
        }

        // Now allocate curr_diag for this diagonal
        int diag_len = i_end - i_start;
        resize_fill_neg_inf(curr_diag, diag_len);

        // We will compute i_start_curr = i_start (the offset for indexing).
        // Similarly we have i_start_prev, i_start_prev2 from the previous diagonals.

        // Shift diagonal buffers:
        prev2_diag.swap(prev_diag); // old 'prev' becomes the new 'prev2'
        prev_diag.swap(curr_diag);  // old 'curr' becomes the new 'prev'
        // Now we must fill curr_diag again (which is empty):
        resize_fill_neg_inf(curr_diag, diag_len);

        // Update the i_start offsets
        i_start_prev2 = i_start_prev;
        i_start_prev  = i_start; 

        // fill the new diagonal
        for(int i = i_start; i < i_end; i++) {
            int j = d - i;
            // boundaries
            if(j < 0 || j > lenB) {
                continue;
            }
            // handle boundary conditions
            if(i == 0 || j == 0) {
                // DP[i][0] = -i, DP[0][j] = -j
                double val = boundary_condition(i, j);
                int idx_curr = i - i_start; 
                curr_diag[idx_curr] = val;
                if(val > max_score) max_score = val;
                explored_count++;
                continue;
            }
            // If not boundary, compute from neighbors
            double cell_val = compute_cell_score(
                i, d,
                prev2_diag,     // diag(d-2)
                prev_diag,      // diag(d-1)
                curr_diag,      // diag(d)
                i_start_prev2,
                i_start_prev,
                i_start
            );
            explored_count++;
            if(cell_val > max_score) {
                max_score = cell_val;
            }
        }

        // Now apply x-drop pruning to the current diagonal
        xdrop_prune(curr_diag);

        // Trim leading/trailing -∞ to shrink i_start/i_end
        trim_i_start(curr_diag, i_start, i_end);
        trim_i_end(curr_diag,   i_start, i_end);

        // If after trimming there's no valid range, we can break
        if(i_end <= i_start) {
            break;
        }
        final_ad = d;
    }
    printf("final ad is %d\n", final_ad);
    // best_score is the maximum cell value found
    return {max_score, explored_count};
}

int main()
{
    // Example: Align random strings
    std::string A = string_gen(10000);
    std::string B = string_gen(10000);

    double x_thresh = 2.0;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto [final_score, visited_count] = xdrop_nw_3diag(A, B, x_thresh);
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_sec = 
        std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    std::cout << "Final alignment score: " << final_score << "\n";
    std::cout << "Cells visited: " << visited_count << "\n";
    std::cout << "Time: " << elapsed_sec << " s\n";
    return 0;
}
