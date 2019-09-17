#include <tmpc/blasfeo/Blasfeo.hpp>

#include <benchmark/benchmark.h>

#include <random>
#include <memory>


#define ADD_BM_GEMM(m, n, p) BENCHMARK_CAPTURE(BM_gemm, m##x##n##x##p##_blasfeo, m, n, p)


namespace tmpc :: benchmark
{
    template <typename MT>
    static void randomize(blasfeo::Matrix<MT>& A)
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(-1.0, 1.0);	

        for (size_t i = 0; i < rows(~A); ++i)
            for (size_t j = 0; j < columns(~A); ++j)
                (~A)(i, j) = dis(gen);
    }


    static void BM_gemm(::benchmark::State& state, size_t m, size_t n, size_t k)
    {
        blasfeo::DynamicMatrix<double> A(k, m), B(n, k), C(m, n), D(m, n);

        randomize(A);
        randomize(B);
        randomize(C);

        /// @brief D <= beta * C + alpha * A^T * B
        // inline void gemm_tn(size_t m, size_t n, size_t k,
        //     double alpha,
        //     blasfeo_dmat const& sA, size_t ai, size_t aj,
        //     blasfeo_dmat const& sB, size_t bi, size_t bj,
        //     double beta,
        //     blasfeo_dmat const& sC, size_t ci, size_t cj,
        //     blasfeo_dmat& sD, size_t di, size_t dj);
        
        for (auto _ : state)
            gemm_tn(m, n, k, 1., A, 0, 0, B, 0, 0, 1., C, 0, 0, D, 0, 0);
    }
    

    ADD_BM_GEMM(2, 2, 2);
    ADD_BM_GEMM(3, 3, 3);
    ADD_BM_GEMM(5, 5, 5);
    ADD_BM_GEMM(10, 10, 10);
    ADD_BM_GEMM(20, 20, 20);
    ADD_BM_GEMM(30, 30, 30);
}
