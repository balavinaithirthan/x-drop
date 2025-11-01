#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

inline int score(char a, char b)
{
    return (a == b) ? 1 : -1;
}

struct AlignmentResult
{
    string aligned_query;
    string aligned_ref;
    int score;
};

// ---------- Normal Correct Backtracking ----------
AlignmentResult backtrack_correct(const vector<vector<int>> &S,
                                  const string &query, const string &ref,
                                  int gap_penalty)
{
    int i = query.size();
    int j = ref.size();

    string A = "", B = "";
    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0 &&
            S[i][j] == S[i - 1][j - 1] + score(query[i - 1], ref[j - 1]))
        {
            A = query[i - 1] + A;
            B = ref[j - 1] + B;
            --i;
            --j;
        }
        else if (i > 0 && S[i][j] == S[i - 1][j] + gap_penalty)
        {
            A = query[i - 1] + A;
            B = '-' + B;
            --i;
        }
        else
        {
            A = '-' + A;
            B = ref[j - 1] + B;
            --j;
        }
    }
    return {A, B, S[query.size()][ref.size()]};
}

// append_function
void append(string &A, char c)
{
    A = c + A; // prepend to front
}

// ---------- Your Alternative Backtracking ----------
AlignmentResult backtrack_alt(const vector<vector<int>> &S,
                              const string &query, const string &ref,
                              int gap_penalty)
{
    int i = query.size();
    int j = ref.size();

    string A = "", B = "";
    std::cout << S[0][0] << "\n";
    while (1)
    {
        if (((i == 0) && (j == 0)))
        {
            break;
        }
        if (((i == i) && (j == 0)))
        {
            if ((S[i][j] == (gap_penalty * i)))
            {
                append(A, query[i - 1]);
                append(B, 0);
                i = (i - 1);
                j = (j + 0);
                continue;
            }
        }
        if (((i == 0) && (j == j)))
        {
            if ((S[i][j] == (gap_penalty * j)))
            {
                append(A, 0);
                append(B, ref[j - 1]);
                i = (i + 0);
                j = (j - 1);
                continue;
            }
        }
        if (((i == i) && (j == j)))
        {
            if ((S[i][j] == (S[i - 1][j - 1] + score(query[i - 1], ref[j - 1]))))
            {
                append(A, query[i - 1]);
                append(B, ref[j - 1]);
                i = (i + -1);
                j = (j + -1);
                continue;
            }
            if ((S[i][j] == (S[i - 1][j] + gap_penalty)))
            {
                append(A, query[i - 1]);
                append(B, '-');
                i = (i - 1);
                j = (j);
                continue;
            }
            if ((S[i][j] == (S[i][j - 1] + gap_penalty)))
            {
                append(A, '-');
                append(B, ref[j - 1]);
                i = (i);
                j = (j - 1);
                continue;
            }
        }
    }
    return {A, B, S[query.size()][ref.size()]};
}

// ---------- Needleman-Wunsch Algorithm ----------
vector<vector<int>> needleman_wunsch(const string &query, const string &ref, int gap_penalty)
{
    int m = query.size(), n = ref.size();
    vector<vector<int>> S(m + 1, vector<int>(n + 1, 0));

    // init
    for (int i = 1; i <= m; i++)
        S[i][0] = i * gap_penalty;
    for (int j = 1; j <= n; j++)
        S[0][j] = j * gap_penalty;

    // fill
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            int diag = S[i - 1][j - 1] + score(query[i - 1], ref[j - 1]);
            int up = S[i - 1][j] + gap_penalty;
            int left = S[i][j - 1] + gap_penalty;
            S[i][j] = max({diag, up, left});
        }
    }

    return S;
}

// ---------- Test Function for Single Example ----------
void test_alignment(const string &query, const string &ref, int gap_penalty, int example_num)
{
    cout << "=== Example " << example_num << " ===\n";
    cout << "Query: " << query << "\n";
    cout << "Ref  : " << ref << "\n";
    cout << "Gap penalty: " << gap_penalty << "\n\n";

    vector<vector<int>> S = needleman_wunsch(query, ref, gap_penalty);

    AlignmentResult correct = backtrack_correct(S, query, ref, gap_penalty);
    AlignmentResult alt = backtrack_alt(S, query, ref, gap_penalty);

    cout << "--- Correct Backtracking ---\n";
    cout << "Aligned Query: " << correct.aligned_query << "\n";
    cout << "Aligned Ref  : " << correct.aligned_ref << "\n";
    cout << "Score: " << correct.score << "\n\n";

    cout << "--- Alternative Backtracking ---\n";
    cout << "Aligned Query: " << alt.aligned_query << "\n";
    cout << "Aligned Ref  : " << alt.aligned_ref << "\n";
    cout << "Score: " << alt.score << "\n\n";

    cout << "Match: " << (correct.aligned_query == alt.aligned_query && correct.aligned_ref == alt.aligned_ref ? "YES" : "NO") << "\n";
    cout << "----------------------------------------\n\n";
}

// ---------- Main with Multiple Examples ----------
int main()
{
    // Example 1: Original test case
    test_alignment("GATTACA", "GCATGCU", -2, 1);

    // Example 2: Simple identical strings
    test_alignment("ATCG", "ATCG", -1, 2);

    // Example 3: Completely different strings
    test_alignment("AAA", "TTT", -1, 3);

    // Example 4: One string is substring of another
    test_alignment("GAT", "GATTACA", -2, 4);

    // Example 5: Empty vs non-empty (edge case)
    test_alignment("", "ABC", -1, 5);

    // Example 6: Longer sequences
    test_alignment("ACGTACGT", "AGCTACGT", -1, 6);

    return 0;
}
