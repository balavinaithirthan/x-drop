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
void aligned_print(std::string query_name, std::string query, std::string ref_name, std::string reference)
{
    // Print the names
    // std::cout << query_name << ": ";

    // Print query sequence with spaces between characters
    for (size_t i = 0; i < query.length(); ++i)
    {
        std::cout << query[i];
        if (i < query.length() - 1)
        {
            std::cout << " ";
        }
    }
    std::cout << std::endl;
    // std::cout << ref_name << ": ";

    // Print reference sequence with spaces between characters
    for (size_t i = 0; i < reference.length(); ++i)
    {
        std::cout << reference[i];
        if (i < reference.length() - 1)
        {
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}
// TODO, have more modularity than just appending to the front
void append(std::vector<char> &A, char c) { A.insert(A.begin(), c); }
int main()
{
    std::string alpha = "ATCG";
    std::vector<int> flat = {
        3, -2, -2, -2, -2, 3, -2, -2, -2, -2, 3, -2, -2, -2, -2, 3};
    SubstitutionMatrix<int> sigma(alpha, flat);
    std::string query_string = "TTTTACCTGAGGGG";
    std::string reference_string = "AAAAACGTAAAAA";
    int M = query_string.size();
    int N = reference_string.size();
    int K = M + N;
    std::vector<size_t> shape = {(size_t)(K + 1), (size_t)(N + 1)};
    Matrix<int> S(shape, 0);
    std::vector<size_t> shapequery = {query_string.size()};
    Matrix<char> query(shapequery, std::vector<char>(query_string.begin(), query_string.end()));
    std::vector<size_t> shapereference = {reference_string.size()};
    Matrix<char> reference(shapereference, std::vector<char>(reference_string.begin(), reference_string.end()));
    std::vector<char> query_aligned;
    std::vector<char> ref_aligned;
    S.print_shape();
    {
        for (int k = 0; k <= K; ++k)
        {
            for (int m = std::max(0, k - M); m <= std::min(N, k); ++m)
            {
                if (((k == 0) && (m == 0)))
                {
                    S(0, 0) = 0;
                    continue;
                }
                if (((k == k - m) && (m == 0)))
                {
                    S((k - m), 0) = (-1 * (k - m));
                    continue;
                }
                if ((k == m) && (m == m))
                {
                    S(m, m) = (-1 * m);
                    continue;
                }
                S(k, m) = variadic_max((S((k - 2), (m - 1)) + sigma(query((k - m - 1)), reference(m - 1))), (S((k - 1), m) + -1), (S((k - 1), (m - 1)) + -1));
                continue;
            }
        }
        for (int k = 0; k <= K; ++k)
        {
            // Print matrix values with proper spacing
            for (int m = std::max(0, k - M); m <= std::min(N, k); ++m)
            // for (size_t m = 0; m < N + 1; ++m)
            {
                std::cout << std::setw(4) << S(k, m);
            }
            std::cout << "\n";
        }
    }
}