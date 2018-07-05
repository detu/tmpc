#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include <benchmark/benchmark.h>

#include <random>
#include <memory>


#define ADD_BM_MTIMES(m, n, p) BENCHMARK_CAPTURE(BM_MTimes, m##x##n##x##p##_blasfeo, m, n, p)


namespace tmpc :: benchmark
{
    static void randomize(size_t m, size_t n, blasfeo_dmat * A)
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(-1.0, 1.0);	

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                BLASFEO_DMATEL(A, i, j) = dis(gen);
    }


    static auto alignedAlloc(size_t bytes)
    {
        return std::unique_ptr<char[], decltype(&std::free)>(
            reinterpret_cast<char *>(std::aligned_alloc(0x1000, bytes)), &std::free);
    }


    static void BM_MTimes(::benchmark::State& state, size_t m, size_t n, size_t p)
    {
        blasfeo_dmat A, B, C;

        auto mem_A = alignedAlloc(blasfeo_memsize_dmat(m, n));
        auto mem_B = alignedAlloc(blasfeo_memsize_dmat(n, p));
        auto mem_C = alignedAlloc(blasfeo_memsize_dmat(m, p));

        blasfeo_create_dmat(m, n, &A, mem_A.get());
        blasfeo_create_dmat(n, p, &B, mem_B.get());
        blasfeo_create_dmat(m, p, &C, mem_C.get());

        randomize(m, n, &A);
        randomize(n, p, &B);

        //blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);
        
        for (auto _ : state)
            // m, n: "Rows and columns of B according to the netlib docs" according to @zanellia
            blasfeo_dtrmm_rlnn(n /*m*/, p /*n*/, 1., &A, 0, 0, &B, 0, 0, &C, 0, 0);
    }

    
    ADD_BM_MTIMES(2, 2, 2);
    ADD_BM_MTIMES(3, 3, 3);
    ADD_BM_MTIMES(5, 5, 5);
    ADD_BM_MTIMES(10, 10, 10);
    ADD_BM_MTIMES(20, 20, 20);
    ADD_BM_MTIMES(30, 30, 30);
}
