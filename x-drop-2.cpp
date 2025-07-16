
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip> // for std::setw

std::string query = "ACGTACGA";    // query sequence
std::string sequence = "ACGTAAAT"; // sequence
const int rows = 9;                // rows
const int cols = 9;                // cols
const int m = rows - 1;            // number of rows in the DP matrix
const int n = cols - 1;            // number of columns in the DP matrix
const int gap = -10;               // gap penalty
const int match = 10;              // match score
const int mismatch = -5;           // mismatch score
const int NEG_INF = -10000;        // negative infinity for initialization
const int INF = 10000;             // positive infinity for initialization
int x_drop = 5;

int two_to_one_dim(int d, int i)
{
  int j = d - i; // j is the second dimension index
  return i * cols + j;
}

int match_fn(char a, char b)
{
  if (a == b)
  {
    return match; // match score
  }
  else
  {
    return mismatch; // mismatch score
  }
}

void print_matrix(const std::vector<int> &S, int m, int n)
{
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      int val = S[i * cols + j];
      if (val == -10000)
      {
        std::cout << std::setw(6) << ".";
      }
      else
      {
        std::cout << std::setw(6) << val;
      }
    }
    std::cout << "\n";
  }
}

int needleman_wunsch()
{
  // g++ -std=c++17 x-drop-2.cpp
  int max_score = 0;
  int var_lo = NEG_INF;
  int var_hi = INF; // var_lo < d - 2i < var_hi
  std::vector<int> S(rows * cols, -1000);
  for (int d = 0; d < ((m + n) + 1); ++d)
  {
    int i_start = std::max(0, d - n);
    // i_start = std::max(i_start, (d - INF) / 2); // we want i_start to always grow
    i_start = std::max(i_start, (d - var_hi) / 2); // we want i_start to always grow

    int i_end = std::min(m, d);
    // i_end = std::min(i_end, (d - NEG_INF) / 2);
    i_end = std::min(i_end, (d - var_lo) / 2);
    // printf("i_start = %d, i_end = %d, d = %d\n", i_start, i_end, d);
    // int i_start = std::max(((var_lo) / 2) + (d / 2), std::max(0, (d - n)));
    // int i_end = std::min((((-var_hi) / 2) + (d / 2)), std::min(m, d));
    // i_start = std::max(0, d - n);
    // i_end = std::min(m, d);
    for (int i = i_start; i <= i_end; ++i)
    {
      // printf("query[%d] = %c, sequence[%d] = %c\n", i, query[i], d - i, sequence[d - i]);
      if (i == 0)
      {
        S[two_to_one_dim(d, i)] = (d - i) * gap;
        if ((S[two_to_one_dim(d, i)] < (max_score - x_drop)))
        {
          S[two_to_one_dim(d, i)] = NEG_INF;
        }
      }
      else if ((d - i) == 0)
      {
        S[two_to_one_dim(d, i)] = i * gap;
        if ((S[two_to_one_dim(d, i)] < (max_score - x_drop)))
        {
          S[two_to_one_dim(d, i)] = NEG_INF;
        }
      }
      else
      {
        S[two_to_one_dim(d, i)] = std::max(
            S[two_to_one_dim(d - 1, i - 1)] + match_fn(query.at(i - 1), sequence.at(d - i - 1)), // match
            std::max(
                S[two_to_one_dim(d - 1, i)] - gap,    // delete
                S[two_to_one_dim(d - 2, i - 1)] - gap // insert
                ));
        // S[ad−1, i−1] + 1 S[ad−1, i] + 1 S[ad−2, i−1] + δ

        if ((S[two_to_one_dim(d, i)] < (max_score - x_drop)))
        {
          S[two_to_one_dim(d, i)] = NEG_INF;
        }
      }
      max_score = std::max(S[two_to_one_dim(d, i)], max_score);
      // S[two_to_one_dim(d, i)] = 5;
    }
    // var_lo < d - 2i < var_hi
    int last_good = i_end;
    for (int i = i_end; i >= i_start; i--) // shrinks e_end
    {
      {
        if ((S[two_to_one_dim(d, i)] != NEG_INF))
        {
          // printf("hi \n");
          last_good = i;
          break;
        }
      } // remember i_end = min(i_end, (d - var_lo) / 2); and that var_lo will be negative
    }
    if (last_good != i_end)
    {
      // printf("last_good = %d, d = %d\n", last_good, d);
      var_lo = std::max(((d - last_good) - last_good), var_lo);
    }
    last_good = i_start;
    for (int i = i_start; i <= i_end; ++i)
    {
      if ((S[two_to_one_dim(d, i)] != NEG_INF))
      {
        last_good = i;
        break;
      }
    }
    if (last_good != i_start)
    {
      printf("last_good = %d, d = %d\n", last_good, d);

      var_hi = std::min(((d - last_good) - last_good), var_hi);
    }
    printf("var_lo = %d, var_hi = %d\n", var_lo, var_hi);
  }
  print_matrix(S, m, n);
  return 0;
}
int main()
{
  int score = needleman_wunsch();
  std::cout << "\nFinal alignment score: " << score << "\n";
  return 0;
}
