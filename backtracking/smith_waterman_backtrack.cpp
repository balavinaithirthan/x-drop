#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

void printMatrix(const vector<vector<int>> &matrix, const string &seq1, const string &seq2)
{
    cout << "   ";
    for (int j = 0; j < seq2.length(); j++)
    {
        cout << "  " << seq2[j];
    }
    cout << endl;

    for (int i = 0; i < matrix.size(); i++)
    {
        if (i == 0)
            cout << " ";
        else
            cout << seq1[i - 1];

        for (int j = 0; j < matrix[i].size(); j++)
        {
            cout << "  " << matrix[i][j];
        }
        cout << endl;
    }
}

int main()
{
    // Input sequences
    string seq1 = "ACGT";
    string seq2 = "AGT";

    cout << "Smith-Waterman Algorithm - Simple Implementation" << endl;
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "Match: +2, Mismatch: -1, Gap: -1" << endl
         << endl;

    // Scoring parameters
    int match = 2;
    int mismatch = -1;
    int gap = -1;

    // Create matrix (m+1) x (n+1)
    int m = seq1.length();
    int n = seq2.length();
    vector<vector<int>> matrix(m + 1, vector<int>(n + 1, 0));

    // Initialize first row and column with zeros (Smith-Waterman property)
    for (int i = 0; i <= m; i++)
    {
        matrix[i][0] = 0;
    }
    for (int j = 0; j <= n; j++)
    {
        matrix[0][j] = 0;
    }

    // Fill the matrix
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // Calculate scores for three possible moves
            int diagonal_score = matrix[i - 1][j - 1];
            if (seq1[i - 1] == seq2[j - 1])
            {
                diagonal_score += match;
            }
            else
            {
                diagonal_score += mismatch;
            }

            int up_score = matrix[i - 1][j] + gap;
            int left_score = matrix[i][j - 1] + gap;

            // Take maximum of the three scores, but not less than 0 (Smith-Waterman)
            matrix[i][j] = max(0, max(diagonal_score, max(up_score, left_score)));
        }
    }

    // Print the filled matrix
    cout << "Filled Matrix:" << endl;
    printMatrix(matrix, seq1, seq2);

    // Find maximum score and its position
    int max_score = 0;
    int max_i = 0, max_j = 0;
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            if (matrix[i][j] > max_score)
            {
                max_score = matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    cout << endl
         << "Maximum score: " << max_score << " at position (" << max_i << ", " << max_j << ")" << endl;

    auto start_i = max_i;
    auto start_j = max_j;
    std::string output_seq1, output_seq2;

    // start position
    // end condition
    // update for each output
    while (true)
    {
        if (start_i == 0 || start_j == 0 || matrix[start_i][start_j] == 0)
            break;

        if (matrix[start_i][start_j] == matrix[start_i - 1][start_j - 1] + (seq1[start_i - 1] == seq2[start_j - 1] ? match : mismatch))
        {
            output_seq1 += seq1[start_i - 1];
            output_seq2 += seq2[start_j - 1];
            start_i--;
            start_j--;
        }
        else if (matrix[start_i][start_j] == matrix[start_i - 1][start_j] + gap)
        {
            output_seq1 += seq1[start_i - 1];
            output_seq2 += '-';
            start_i--;
        }
        else
        {
            output_seq1 += '-';
            output_seq2 += seq2[start_j - 1];
            start_j--;
        }
    }

    return 0;
}
