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
#include <chrono>

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
        {
            std::cerr << "Unknown symbol: " << a << " or " << b << std::endl;
            throw std::out_of_range("unknown symbol");
        }
        return data_[it1->second * N_ + it2->second];
    }
};

template <typename T>
class Matrix
{
public:
    // Constructor with fill value (existing functionality)
    Matrix(const std::vector<size_t> &shape, const T &fill_value)
        : shape_(shape)
    {
        size_t total_size = std::accumulate(shape.begin(), shape.end(), 1UL,
                                            std::multiplies<size_t>());
        data_ = std::vector<T>(total_size, fill_value);

        // Compute strides
        strides_.resize(shape.size());
        strides_.back() = 1;
        for (int i = static_cast<int>(shape.size()) - 2; i >= 0; --i)
        {
            strides_[i] = strides_[i + 1] * shape[i + 1];
        }
    }

    // Constructor with optional data vector
    Matrix(const std::vector<size_t> &shape, const std::vector<T> &data)
        : shape_(shape)
    {
        size_t total_size = std::accumulate(shape.begin(), shape.end(), 1UL,
                                            std::multiplies<size_t>());

        if (data.empty())
        {
            // If no data provided, initialize with default value
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

        // Compute strides
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
            if (i != shape_.size() - 1)
                std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }

    std::vector<size_t> shape() const
    {
        return shape_;
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

Matrix<int> original_align_sequences(const std::vector<char> &query_seq, const std::vector<char> &ref_seq)
{
    std::string alpha = "ATCG";
    std::vector<int> flat = {
        2, -1, -1, -1,
        -1, 2, -1, -1,
        -1, -1, 2, -1,
        -1, -1, -1, 2};
    SubstitutionMatrix<int> sigma(alpha, flat);

    size_t K = query_seq.size() + ref_seq.size(); // max combined length
    size_t M = ref_seq.size();

    Matrix<int> S({K + 1, M + 1}, 0); // DP matrix

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t k = 0; k <= K; ++k)
    {
        for (size_t m = std::max<int>(0, k - query_seq.size()); m <= std::min(k, M); ++m)
        {
            if (k == 0 && m == 0)
            {
                S(0, 0) = 0;
            }
            else if (m == 0 && k - m == k)
            {
                S(k, 0) = -2 * k;
            }
            else if (k == m)
            {
                S(m, m) = -2 * m;
            }
            else
            {
                S(k, m) = variadic_max(S(k - 2, m - 1) + sigma(query_seq[k - m - 1], ref_seq[m - 1]),
                                       S(k - 1, m) - 2,
                                       S(k - 1, m - 1) - 2);
            }
        }
    }

    // End timing and print result
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Original algorithm execution time: " << duration.count() << " microseconds" << std::endl;

    return S;
}

Matrix<int> faster_align_sequences(const std::vector<char> &query_seq, const std::vector<char> &ref_seq)
{
    std::string alpha = "ATCG";
    std::vector<int> flat = {
        2, -1, -1, -1,
        -1, 2, -1, -1,
        -1, -1, 2, -1,
        -1, -1, -1, 2};
    SubstitutionMatrix<int> sigma(alpha, flat);

    size_t K = query_seq.size() + ref_seq.size(); // max combined length
    size_t M = ref_seq.size();

    Matrix<int> S({K + 1, M + 1}, 0); // DP matrix

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    S(0, 0) = 0; // Base case
    for (size_t k = 1; k <= K; ++k)
    {
        S(k, 0) = -2 * k; // First column
    }
    for (size_t m = 1; m <= M; ++m)
    {
        S(m, m) = -2 * m; // First row
    }

    for (size_t k = 2; k <= K; ++k)
    {
        for (size_t m = std::max<int>(1, k - query_seq.size()); m <= std::min(k - 1, M); ++m)
        {
            S(k, m) = variadic_max(S(k - 2, m - 1) + sigma(query_seq[k - m - 1], ref_seq[m - 1]),
                                   S(k - 1, m) - 2,
                                   S(k - 1, m - 1) - 2);
        }
    }

    // End timing and print result
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Faster algorithm execution time: " << duration.count() << " microseconds" << std::endl;

    return S;
}

// main
int main()
{
    std::vector<char> query_seq(3000, 'C'); // 200 As
    std::vector<char> ref_seq(1000, 'A');   // 300 As

    Matrix<int> alignment = original_align_sequences(query_seq, ref_seq);
    alignment.print_shape();
    std::cout << "Original Alignment:" << std::endl;
    // print S
    // for (size_t i = 0; i < alignment.shape()[0]; ++i)
    // {
    //     for (size_t j = 0; j < alignment.shape()[1]; ++j)
    //     {
    //         std::cout << alignment(i, j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    Matrix<int> faster_alignment = faster_align_sequences(query_seq, ref_seq);
    faster_alignment.print_shape();
    std::cout << "Faster Alignment:" << std::endl;
    // print S
    // for (size_t i = 0; i < faster_alignment.shape()[0]; ++i)
    // {
    //     for (size_t j = 0; j < faster_alignment.shape()[1]; ++j)
    //     {
    //         std::cout << faster_alignment(i, j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // percent improvement
    return 0;
}

/*
Conclusion

- all mismatches 30% slower
- all matches is the same

*/