#include <algorithm>
#include <fstream>
#include <iomanip> // for std::setw
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

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

int to1D(int i, int j, int num_cols) { return i * num_cols + j; }

/*
 * Initialize an N x N matrix (where N = len(A) + ping) with 'fill_val' in its
 * interior, but sets the base-case boundaries: mat[i][0] = -i, mat[0][j] = -j.
 *
 * NOTE: In the original Python code, for demonstration, the matrix is sized
 * by A’s length alone, ignoring B’s length. This works for the example if A
 * and B are the same length. Adjust for the general case if needed.
 */
std::vector<int> initializeMatrix(const std::string &A, const std::string &B,
                                  int fill_val, int &num_rows, int &num_cols) {
  int N = static_cast<int>(A.size());
  num_rows = N + ping;
  num_cols = N + ping;
  std::vector<int> mat(num_rows * num_cols, 0);

  // Fill interior (excluding the top row & left col) with fill_val
  for (int i = 1; i < num_rows; i++) {
    for (int j = 1; j < num_cols; j++) {
      mat[to1D(i, j, num_cols)] = fill_val;
    }
  }

  // Base case
  for (int i = 1; i < static_cast<int>(A.size()) + ping; i++) {
    mat[to1D(i, 0, num_cols)] = -i;
  }
  for (int j = 1; j < static_cast<int>(B.size()) + ping; j++) {
    if (j < num_cols)
      mat[to1D(0, j, num_cols)] = -j;
  }
  return mat;
}

/**
 * Print a matrix for debugging.
 */
void printMatrix(const std::vector<int> &mat, int num_rows, int num_cols) {
  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      std::cout << std::setw(6) << mat[to1D(i, j, num_cols)] << " ";
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
int score(const std::vector<int> &mat, const std::string &A,
          const std::string &B, int i, int j, int num_cols) {
  int diag = mat[to1D(i - 1, j - 1, num_cols)] + (A[i - 1] == B[j - 1] ? 1 : 0);
  int left = mat[to1D(i, j - 1, num_cols)] - 2;
  int up = mat[to1D(i - 1, j, num_cols)] - 2;
  return std::max({diag, left, up});
}

/**
 * get_i_start() from the Python code.
 *   We move i upward from i_start while mat[i][j] == -∞, then return the new i
 *   and a diagonal offset d_lo.
 */
std::pair<int, int> get_i_start(std::vector<int> &mat,
                                std::vector<int> &explored, int d, int i_start,
                                int i_end, int num_cols) {
  int i = i_start;
  int d_lo = -std::numeric_limits<int>::infinity();
  while (i < i_end) {
    int j = d - i;
    if (j < 0 || j >= num_cols)
      break;
    if (mat[to1D(i, j, num_cols)] != -std::numeric_limits<int>::infinity()) {
      // We found a valid cell
      break;
    }
    // Mark explored
    // explored[to1D(i, j, num_cols)] = 1;
    i++;
  }
  // If the original i_start cell is not -∞, set d_lo to -∞
  // else set something else. (Mirroring the Python code logic.)
  if (i_start < i_end) {
    int j_start = d - i_start;
    if (0 <= j_start && j_start < num_cols) {
      if (mat[to1D(i_start, j_start, num_cols)] !=
          -std::numeric_limits<int>::infinity()) {
        d_lo = -std::numeric_limits<int>::infinity();
      } else {
        int j_temp = d - i;
        d_lo = 2 * i - d;
      }
    }
  }
  return std::make_pair(i, d_lo);
}

/**
 * get_i_end() from Python code. We move downward from i_end-1
 * while mat[i][j] == -∞, then return new i_end and d_hi.
 */
std::pair<int, int> get_i_end(std::vector<int> &mat, std::vector<int> &explored,
                              int d, int i_start, int i_end, int num_cols) {
  int i = i_end - 1; // we walk downward from i_end-1
  int d_hi = std::numeric_limits<int>::infinity();
  while (i >= i_start) {
    int j = d - i;
    if (j < 0 || j >= num_cols)
      break;
    if (mat[to1D(i, j, num_cols)] != -std::numeric_limits<int>::infinity()) {
      // Found a valid cell
      break;
    }
    // explored[to1D(i, j, num_cols)] = 1;
    i--;
  }
  // If mat[i_end-1][?] != -∞, set d_hi = ∞ else something else
  int i_chk = i_end - 1;
  int j_chk = d - i_chk;
  if (i_chk >= 0 && i_chk < num_cols && j_chk >= 0 && j_chk < num_cols) {
    if (mat[to1D(i_chk, j_chk, num_cols)] !=
        -std::numeric_limits<int>::infinity()) {
      d_hi = std::numeric_limits<int>::infinity();
    } else {
      int j_temp = d - i;
      d_hi = 2 * i - d + 2;
    }
  }
  return std::make_pair(i + 1, d_hi);
}

/**
 * Mark cells whose score < (max_score - X_drop) as -∞ and also mark 'explored'.
 */
void markX(std::vector<int> &mat, std::vector<int> &explored, int d,
           int i_start, int i_end, int X_drop, int max_score, int num_cols) {
  for (int i = i_start; i < i_end; i++) {
    int j = d - i;
    if (j < 0 || j >= num_cols)
      continue;
    int s = mat[to1D(i, j, num_cols)];
    // explored[to1D(i, j, num_cols)] = 1;
    if (s < (max_score - X_drop)) {
      mat[to1D(i, j, num_cols)] = -std::numeric_limits<int>::infinity();
    }
  }
}

/**
 * modify_i_start: from python code
 *   i_on_diag = (d_lo + ad) // 2
 *   i_start = max(i_start, i_on_diag)
 */
int modify_i_start(int i_start, int d_lo, int ad) {
  if (d_lo == -std::numeric_limits<int>::infinity()) {
    // same as python: if no real offset, just keep i_start
    return i_start;
  }
  // The python code does (d_lo + ad) // 2 for i_on_diag
  // We'll do integer division in C++ as well:
  int i_on_diag = (d_lo + ad) / 2;
  return std::max(i_start, i_on_diag);
}

/**
 * modify_i_end: from python code
 *   i_on_diag = (d_hi + ad) // 2
 *   i_end = min(i_end, i_on_diag)
 */
int modify_i_end(int i_end, int d_hi, int ad) {
  if (d_hi == std::numeric_limits<int>::infinity()) {
    return i_end;
  }
  int i_on_diag = (d_hi + ad) / 2;
  return std::min(i_end, i_on_diag);
}

/**
 * The main x-drop NW function, replicating the Python's nw4(A, B, x_thresh=2).
 * Returns the final matrix with computed scores and the explored matrix.
 *
 * The logic is specialized for the example, but follows the Python code
 * anti-diagonal iteration.
 */
std::pair<std::vector<int>, std::vector<int>> xdrop_nw(const std::string &A,
                                                       const std::string &B,
                                                       int x_thresh,
                                                       std::ofstream &outfile) {
  int num_rows, num_cols;
  std::vector<int> mat = initializeMatrix(A, B, 1, num_rows, num_cols);
  std::vector<int> explored;
  // initializeMatrix(A, B, 0, num_rows, num_cols);

  auto start_time = std::chrono::high_resolution_clock::now();
  int m = static_cast<int>(A.size()) + 1; // rows
  int n = static_cast<int>(B.size()) + 1; // cols

  int i_start = 1;
  int i_end = 1;
  int max_score = -std::numeric_limits<int>::infinity();

  int lower_diag = -std::numeric_limits<int>::infinity();
  int upper_diag = std::numeric_limits<int>::infinity();
  int final_ad = 0;
  // First pass: ad in [2, m)
  for (int ad = 2; ad < m; ad++) {
    i_end++;
    for (int i = i_start; i < i_end; i++) {
      int j = ad - i;
      if (j < 0 || j >= num_cols) {
        continue;
      }
      // compute the score
      int s = score(mat, A, B, i, j, num_cols);
      mat[to1D(i, j, num_cols)] = s;
      // explored[to1D(i, j, num_cols)] = 1;
      if (s > max_score) {
        max_score = s;
      }
    }
    // Mark X-drop
    markX(mat, explored, ad, i_start, i_end, x_thresh, max_score, num_cols);

    // Adjust i_start, i_end
    auto start_pair = get_i_start(mat, explored, ad, i_start, i_end, num_cols);
    i_start = start_pair.first;
    int d_lo = start_pair.second;

    auto end_pair = get_i_end(mat, explored, ad, i_start, i_end, num_cols);
    i_end = end_pair.first;
    int d_hi = end_pair.second;

    lower_diag = std::max(lower_diag, d_lo);
    upper_diag = std::min(upper_diag, d_hi);

    i_start = modify_i_start(i_start, lower_diag, ad);
    i_end = modify_i_end(i_end, upper_diag, ad);
  }

  // Second pass: ad in [m, m + n - 1)
  // This covers anti-diagonals that extend beyond the first dimension's corner.
  for (int ad = m; ad < m + n - 1; ad++) {
    int expected_i_start = ad - m + 1;
    if (expected_i_start > i_start) {
      i_start = expected_i_start;
    }
    int expected_i_end = m;
    if ((i_end + 1) < expected_i_end) {
      i_end = i_end + 1;
    } else {
      i_end = expected_i_end;
    }
    for (int i = i_start; i < i_end; i++) {
      int j = ad - i;
      if (j < 0 || j >= num_cols || i < 0 || i >= num_rows) {
        continue;
      }
      int s = score(mat, A, B, i, j, num_cols);
      mat[to1D(i, j, num_cols)] = s;
      // explored[to1D(i, j, num_cols)] = 1;
      if (s > max_score) {
        max_score = s;
        // }
      }

      if (i_start == i_end) {
        // no more space to move -> break
        break;
      }

      // Mark X-drop
      markX(mat, explored, ad, i_start, i_end, x_thresh, max_score, num_cols);

      // Adjust i_start, i_end
      auto start_pair =
          get_i_start(mat, explored, ad, i_start, i_end, num_cols);
      i_start = start_pair.first;
      int d_lo = start_pair.second;

      auto end_pair = get_i_end(mat, explored, ad, i_start, i_end, num_cols);
      i_end = end_pair.first;
      int d_hi = end_pair.second;

      lower_diag = std::max(lower_diag, d_lo);
      upper_diag = std::min(upper_diag, d_hi);

      i_start = modify_i_start(i_start, lower_diag, ad);
      i_end = modify_i_end(i_end, upper_diag, ad);
      final_ad = ad;
    }
  }
  printf("final antidiagonal explored is %d\n", final_ad);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  // outfile << "Elapsed time size " << num_cols - 1 << ", " << elapsed.count()
  outfile << elapsed.count() << " \n";
  return std::make_pair(mat, explored);
}

std::string string_gen(const int length) {
  // The alphabet of possible characters
  static const char alphabet[] = {'A', 'G', 'C', 'T'};

  // Set up random number generation
  std::random_device rd;
  std::mt19937 gen(rd()); // Mersenne Twister
  std::uniform_int_distribution<> dist(0, 3);

  // Build a random string
  std::string randomDNA;
  randomDNA.reserve(length); // reserve for efficiency

  for (int i = 0; i < length; ++i) {
    randomDNA.push_back(alphabet[dist(gen)]);
  }
  return randomDNA;
}

void runExperiments(const std::string &filename) {
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  for (int sequenceLength = 100; sequenceLength <= 12000;
       sequenceLength += 100) {
    std::string A = string_gen(sequenceLength);
    std::string B = string_gen(sequenceLength);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto result = xdrop_nw(A, B, 2, outfile);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // outfile << "Sequence Length: " << sequenceLength
    //         << ", Elapsed time: " << elapsed.count() << " s\n";
  }

  outfile.close();
}

int main() {
  // std::string A = string_gen(10000);
  // std::string B = string_gen(10000);

  // auto start_time = std::chrono::high_resolution_clock::now();
  // auto result = xdrop_nw(A, B, 2);
  // auto end_time = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = end_time - start_time;
  // std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  // std::vector<int> mat = result.first;
  // std::vector<int> explored = result.second;

  // int num_rows = static_cast<int>(A.size()) + 1;
  // int num_cols = static_cast<int>(B.size()) + 1;

  // // Print out the final matrix
  // std::cout << "Final x-drop NW matrix:\n";
  // std::cout << mat[to1D(A.size(), B.size(), num_cols)] << std::endl;

  // // Count the number of explored cells where the value is 1
  // int explored_count = 0;
  // for (const auto &cell : explored) {
  //   if (cell == 1) {
  //     explored_count++;
  //   }
  // }

  // std::cout << "Number of explored cells: " << explored_count << std::endl;
  // std::cout << "Totoal cells in the matrix: " << A.size() * B.size()
  //           << std::endl;

  runExperiments("xdrop_runtime_results.txt");

  return 0;
}
