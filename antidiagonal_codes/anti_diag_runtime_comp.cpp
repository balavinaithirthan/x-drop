#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

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
        size_t total_size = std::accumulate(shape.begin(), shape.end(), 1UL,
                                            std::multiplies<size_t>());
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
        size_t total_size = std::accumulate(shape.begin(), shape.end(), 1UL,
                                            std::multiplies<size_t>());

        if (data.empty())
        {
            data_ = std::vector<T>(total_size, T{});
        }
        else
        {
            if (data.size() != total_size)
            {
                throw std::invalid_argument(
                    "Data size does not match matrix dimensions");
            }
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

    void fill(const T &val)
    {
        std::fill(data_.begin(), data_.end(), val);
    }

private:
    std::vector<size_t> shape_;
    std::vector<size_t> strides_;
    std::vector<T> data_;

    size_t flat_index(const std::vector<size_t> &indices) const
    {
        if (indices.size() != shape_.size())
            throw std::invalid_argument(
                "Number of indices does not match matrix dimension");

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
T variadic_max(T first, Args... args)
{
    return (std::max)({first, args...});
}

template <typename T, typename... Args>
T variadic_min(T first, Args... args)
{
    return (std::min)({first, args...});
}

void aligned_print(std::string /*query_name*/, std::string query,
                   std::string /*ref_name*/, std::string reference)
{
    for (size_t i = 0; i < query.length(); ++i)
    {
        std::cout << query[i] << (i + 1 < query.length() ? " " : "");
    }
    std::cout << "\n";
    for (size_t i = 0; i < reference.length(); ++i)
    {
        std::cout << reference[i] << (i + 1 < reference.length() ? " " : "");
    }
    std::cout << "\n";
}

void append(std::vector<char> &A, char c) { A.insert(A.begin(), c); }

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::string alpha = "ATCG";
    std::vector<int> flat = {
        3, -2, -2, -2,
        -2, 3, -2, -2,
        -2, -2, 3, -2,
        -2, -2, -2, 3};
    SubstitutionMatrix<int> sigma(alpha, flat);

    // std::string query_string = "TTTTACCTGAGGGG";
    // std::string reference_string = "AAAAACGTAAAAA";
    std::string query_string = std::string(4000, 'A');
    std::string reference_string = std::string(4000, 'A');

    int M = static_cast<int>(query_string.size());
    int N = static_cast<int>(reference_string.size());
    int K = M + N;

    Matrix<char> query({query_string.size()}, std::vector<char>(query_string.begin(), query_string.end()));
    Matrix<char> reference({reference_string.size()}, std::vector<char>(reference_string.begin(), reference_string.end()));

    // -----------------------------
    // 1) Antidiagonal (k, m) version
    // -----------------------------
    Matrix<int> S_km({static_cast<size_t>(K + 1), static_cast<size_t>(N + 1)}, 0);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k <= K; ++k)
    {
        int m_start = std::max(0, k - M);
        int m_end = std::min(N, k);

// Parallelize across the m's of the current antidiagonal.
// Dependencies are only to prior diagonals (k-1 or k-2), so this is safe.
#pragma omp parallel for schedule(static)
        for (int m = m_start; m <= m_end; ++m)
        {
            int i = k - m;
            int j = m;

            if (i == 0 && j == 0)
            {
                S_km(0, 0) = 0;
                continue;
            }
            if (j == 0)
            {
                S_km(i, 0) = -1 * i;
                continue;
            }
            if (i == 0)
            {
                S_km(0, j) = -1 * j;
                continue;
            }

            // Recurrence in (k, m) space:
            // i = k - m, j = m
            // S(i,j) = max(
            //   S(i-1, j-1) + sigma(q[i-1], t[j-1]),
            //   S(i-1, j)   - 1,
            //   S(i,   j-1) - 1)
            // Mapping:
            //   (i-1,j-1) -> (k-2, m-1)
            //   (i-1,j)   -> (k-1, m)
            //   (i,  j-1) -> (k-1, m-1)
            int match = S_km(k - 2, m - 1) + sigma(query(i - 1), reference(j - 1));
            int del = S_km(k - 1, m) + -1;
            int insert = S_km(k - 1, m - 1) + -1;
            S_km(k, m) = variadic_max(match, del, insert);
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms_km = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // -----------------------------
    // 2) Plain (i, j) version
    // -----------------------------
    Matrix<int> S_ij({query_string.size() + 1, reference_string.size() + 1}, 0);

    auto t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < static_cast<int>(query_string.size()) + 1; ++i)
    {
        for (int j = 0; j < static_cast<int>(reference_string.size()) + 1; ++j)
        {
            if (i == 0 && j == 0)
            {
                S_ij(0, 0) = 0;
                continue;
            }
            if (j == 0)
            {
                S_ij(i, 0) = (-1 * i);
                continue;
            }
            if (i == 0)
            {
                S_ij(0, j) = (-1 * j);
                continue;
            }
            S_ij(i, j) = variadic_max(
                S_ij(i - 1, j - 1) + sigma(query(i - 1), reference(j - 1)),
                S_ij(i - 1, j) + -1,
                S_ij(i, j - 1) + -1);
        }
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    double ms_ij = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // -----------------------------
    // Report
    // -----------------------------
    std::cout << std::fixed << std::setprecision(3);
#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n";
#else
    std::cout << "OpenMP threads: (not compiled with -fopenmp)\n";
#endif
    std::cout << "Antidiagonal (k,m) time: " << ms_km << " ms\n";
    std::cout << "Row/Col (i,j) time    : " << ms_ij << " ms\n";

    std::cout << "Speedup (i,j)/(k,m)   : " << (ms_ij / ms_km) << "x\n";

    // Optional: show last cell scores to ensure both compute same end state
    // for global alignment.
    // std::cout << "S_ij(M,N) = " << S_ij(query_string.size(), reference_string.size()) << "\n";
    // // Map (M,N) -> (k= M+N, m=N)
    // std::cout << "S_km(K,N) = " << S_km(M + N, N) << "\n";

    return 0;
    // clang++ -O3 -march=native -fopenmp antidiagonal_codes/anti_diag_runtime_comp.cpp -o dp_compare
