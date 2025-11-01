#include <arm_neon.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <cstdlib>

// Scalar version
void add_scalar(const float *a, const float *b, float *c, int n)
{
    for (int i = 0; i < n; ++i)
        c[i] = a[i] + b[i];
}

// NEON vectorized version
void add_neon(const float *a, const float *b, float *c, int n)
{
    int i = 0;
    for (; i <= n - 4; i += 4)
    {
        float32x4_t va = vld1q_f32(a + i);
        float32x4_t vb = vld1q_f32(b + i);
        float32x4_t vc = vaddq_f32(va, vb);
        vst1q_f32(c + i, vc);
    }
    for (; i < n; ++i)
        c[i] = a[i] + b[i];
}

int main()
{
    const int N = 1 << 20; // ~1 million elements
    std::vector<float> a(N), b(N), c(N);
    for (int i = 0; i < N; ++i)
    {
        a[i] = std::rand() / float(RAND_MAX);
        b[i] = std::rand() / float(RAND_MAX);
    }

    constexpr int REPS = 100;
    using clock = std::chrono::high_resolution_clock;

    auto measure = [&](auto func, const char *name)
    {
        auto start = clock::now();
        for (int r = 0; r < REPS; ++r)
            func(a.data(), b.data(), c.data(), N);
        auto end = clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        std::cout << name << ": " << ms << " ms (" << (ms / REPS) << " ms per run)\n";
    };

    measure(add_scalar, "Scalar");
    measure(add_neon, "NEON");

    return 0;
}

// clang++ -O3 simd_arm.cpp -o vector_add
