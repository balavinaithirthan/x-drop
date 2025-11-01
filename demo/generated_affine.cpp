#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

template <typename ScoreT = int>
class SubstitutionMatrix
{
    std::unordered_map<char, size_t> index_;
    std::vector<ScoreT> data_;
    size_t N_;

public:
    SubstitutionMatrix(const std::string &alphabet,
                       const std::vector<ScoreT> &matrix_flat)
        : N_(alphabet.size()), data_(matrix_flat)
    {
        assert(matrix_flat.size() == N_ * N_);
        for (size_t i = 0; i < N_; ++i)
        {
            char c = alphabet[i];
            if (index_.count(c))
                throw std::invalid_argument("duplicate alphabet char");
            index_[c] = i;
        }
    }

    ScoreT operator()(char a, char b) const
    {
        auto it1 = index_.find(a);
        auto it2 = index_.find(b);
        if (it1 == index_.end() || it2 == index_.end())
            throw std::out_of_range("unknown symbol");
        return data_[it1->second * N_ + it2->second];
    }
};

template <typename T>
class Matrix
{
public:
    Matrix(const std::vector<size_t> &shape, const T &fill_value)
        : shape_(shape)
    {
        size_t total_size =
            std::accumulate(shape.begin(), shape.end(), 1UL, std::multiplies<size_t>());
        data_ = std::vector<T>(total_size, fill_value);

        strides_.resize(shape.size());
        strides_.back() = 1;
        for (int i = static_cast<int>(shape.size()) - 2; i >= 0; --i)
        {
            strides_[i] = strides_[i + 1] * shape[i + 1];
        }
    }

    Matrix(const std::vector<size_t> &shape, const std::vector<T> &data)
        : shape_(shape)
    {
        size_t total_size =
            std::accumulate(shape.begin(), shape.end(), 1UL, std::multiplies<size_t>());
        if (data.empty())
            data_ = std::vector<T>(total_size, T{});
        else
        {
            if (data.size() != total_size)
                throw std::invalid_argument("Data size does not match matrix dimensions");
            data_ = data;
        }
        strides_.resize(shape.size());
        strides_.back() = 1;
        for (int i = static_cast<int>(shape.size()) - 2; i >= 0; --i)
        {
            strides_[i] = strides_[i + 1] * shape[i + 1];
        }
    }

    template <typename... Args>
    T &operator()(Args... args)
    {
        std::vector<size_t> indices{static_cast<size_t>(args)...};
        return data_.at(flat_index(indices));
    }
    template <typename... Args>
    const T &operator()(Args... args) const
    {
        std::vector<size_t> indices{static_cast<size_t>(args)...};
        return data_.at(flat_index(indices));
    }

    void print_shape() const
    {
        std::cout << "Shape: (";
        for (size_t i = 0; i < shape_.size(); ++i)
        {
            std::cout << shape_[i];
            if (i + 1 < shape_.size())
                std::cout << ", ";
        }
        std::cout << ")\n";
    }

private:
    std::vector<size_t> shape_;
    std::vector<size_t> strides_;
    std::vector<T> data_;

    size_t flat_index(const std::vector<size_t> &indices) const
    {
        if (indices.size() != shape_.size())
            throw std::invalid_argument("Number of indices does not match matrix dimension");
        size_t idx = 0;
        for (size_t i = 0; i < shape_.size(); ++i)
        {
            if (indices[i] >= shape_[i])
                throw std::out_of_range("Index out of bounds");
            idx += indices[i] * strides_[i];
        }
        return idx;
    }
};

template <typename T, typename... Args>
T variadic_max(T first, Args... args) { return (std::max)({first, args...}); }

template <typename T, typename... Args>
T variadic_min(T first, Args... args) { return (std::min)({first, args...}); }

static inline void append(std::vector<char> &A, char c) { A.insert(A.begin(), c); }

static inline void print_spaced(const std::string &s)
{
    for (size_t i = 0; i < s.size(); ++i)
    {
        std::cout << s[i];
        if (i + 1 < s.size())
            std::cout << ' ';
    }
    std::cout << '\n';
}
static inline void aligned_print(const std::string &q_name,
                                 const std::string &q,
                                 const std::string &r_name,
                                 const std::string &r)
{
    print_spaced(q);
    print_spaced(r);
}

int main()
{
    const int gap_open = -2;
    const int gap_extend = -1;
    const int NEG_INF = std::numeric_limits<int>::min() / 4;

    std::string alpha = "ATCG";
    std::vector<int> flat = {
        2, -1, -1, -1,
        -1, 2, -1, -1,
        -1, -1, 2, -1,
        -1, -1, -1, 2};
    SubstitutionMatrix<int> sigma(alpha, flat);

    // Strings for input
    std::string query_string = "TTTTACCTGAGGGG";
    std::string reference_string = "AAAAACGTAAAAA";
    const size_t m = query_string.size();
    const size_t n = reference_string.size();

    // Wrap for indexed access
    Matrix<char> query({m}, std::vector<char>(query_string.begin(), query_string.end()));
    Matrix<char> reference({n}, std::vector<char>(reference_string.begin(), reference_string.end()));

    // DP tables
    Matrix<int> I({m + 1, n + 1}, NEG_INF);
    Matrix<int> D({m + 1, n + 1}, NEG_INF);
    Matrix<int> S({m + 1, n + 1}, NEG_INF);

    // Single DP loop with base cases embedded
    for (size_t i = 0; i <= m; ++i)
    {
        for (size_t j = 0; j <= n; ++j)
        {

            // Corner
            if (i == 0 && j == 0)
            {
                S(0, 0) = 0;
                I(0, 0) = NEG_INF;
                D(0, 0) = NEG_INF;
                continue;
            }

            // First column: only D and S are valid
            if (j == 0)
            {
                // I(i,0) stays invalid
                I(i, 0) = NEG_INF;
                // Recurrence produces correct linear gap: open at i==1, extend otherwise
                int from_open = (i > 0) ? S(i - 1, 0) + gap_open + gap_extend : NEG_INF;
                int from_ext = (i > 0) ? D(i - 1, 0) + gap_extend : NEG_INF;
                D(i, 0) = (std::max)(from_open, from_ext);
                S(i, 0) = D(i, 0);
                continue;
            }

            // First row: only I and S are valid
            if (i == 0)
            {
                // D(0,j) stays invalid
                D(0, j) = NEG_INF;
                int from_open = (j > 0) ? S(0, j - 1) + gap_open + gap_extend : NEG_INF;
                int from_ext = (j > 0) ? I(0, j - 1) + gap_extend : NEG_INF;
                I(0, j) = (std::max)(from_open, from_ext);
                S(0, j) = I(0, j);
                continue;
            }

            // General case
            {
                // I: from left
                int i_from_open = S(i, j - 1) + gap_open + gap_extend;
                int i_from_extend = I(i, j - 1) + gap_extend;
                I(i, j) = (std::max)(i_from_open, i_from_extend);

                // D: from up
                int d_from_open = S(i - 1, j) + gap_open + gap_extend;
                int d_from_extend = D(i - 1, j) + gap_extend;
                D(i, j) = (std::max)(d_from_open, d_from_extend);

                // S: diag or I/D
                int diag = S(i - 1, j - 1) + sigma(query(i - 1), reference(j - 1));
                S(i, j) = variadic_max(diag, I(i, j), D(i, j));
            }
        }
    }

    // Print S
    for (size_t i = 0; i <= m; ++i)
    {
        for (size_t j = 0; j <= n; ++j)
        {
            if (S(i, j) <= NEG_INF / 2)
                std::cout << std::setw(5) << "-inf";
            else
                std::cout << std::setw(5) << S(i, j);
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    S.print_shape();

    // Backtracking
    enum State
    {
        ST_S,
        ST_I,
        ST_D
    };
    size_t bi = m, bj = n;
    int s_end = S(m, n), i_end = I(m, n), d_end = D(m, n);
    State st;
    if (s_end >= i_end && s_end >= d_end)
        st = ST_S;
    else if (i_end >= d_end)
        st = ST_I;
    else
        st = ST_D;

    std::vector<char> q_aln, r_aln;
    while (bi > 0 || bj > 0)
    {
        switch (st)
        {
        case ST_S:
        {
            if (bi > 0 && bj > 0 &&
                S(bi, bj) == S(bi - 1, bj - 1) + sigma(query(bi - 1), reference(bj - 1)))
            {
                append(q_aln, query(bi - 1));
                append(r_aln, reference(bj - 1));
                --bi;
                --bj;
                break;
            }
            if (S(bi, bj) == I(bi, bj))
            {
                st = ST_I;
                break;
            }
            if (S(bi, bj) == D(bi, bj))
            {
                st = ST_D;
                break;
            }
            throw std::runtime_error("Backtrack in S failed");
        }
        case ST_I:
        {
            if (bj > 0 && I(bi, bj) == S(bi, bj - 1) + gap_open + gap_extend)
            {
                append(q_aln, '-');
                append(r_aln, reference(bj - 1));
                --bj;
                st = ST_S;
                break;
            }
            if (bj > 0 && I(bi, bj) == I(bi, bj - 1) + gap_extend)
            {
                append(q_aln, '-');
                append(r_aln, reference(bj - 1));
                --bj;
                break;
            }
            throw std::runtime_error("Backtrack in I failed");
        }
        case ST_D:
        {
            if (bi > 0 && D(bi, bj) == S(bi - 1, bj) + gap_open + gap_extend)
            {
                append(q_aln, query(bi - 1));
                append(r_aln, '-');
                --bi;
                st = ST_S;
                break;
            }
            if (bi > 0 && D(bi, bj) == D(bi - 1, bj) + gap_extend)
            {
                append(q_aln, query(bi - 1));
                append(r_aln, '-');
                --bi;
                break;
            }
            throw std::runtime_error("Backtrack in D failed");
        }
        }
    }

    // Print original and aligned
    std::cout << "\nOriginal:\n";
    aligned_print("query", query_string,
                  "reference", reference_string);
    std::cout << "\nAligned:\n";
    std::string q_out(q_aln.begin(), q_aln.end());
    std::string r_out(r_aln.begin(), r_aln.end());
    aligned_print("query_aln", q_out, "ref_aln", r_out);

    return 0;
}
