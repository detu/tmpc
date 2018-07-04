#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include <benchmark/benchmark.h>

#include <vector>
#include <random>


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


    static void BM_MTimes(::benchmark::State& state, size_t m, size_t n, size_t p)
    {
        blasfeo_dmat A, B, C;

        std::vector<char> mem_A(::blasfeo_memsize_dmat(m, n));
        std::vector<char> mem_B(blasfeo_memsize_dmat(n, p));
        std::vector<char> mem_C(blasfeo_memsize_dmat(m, p));

        blasfeo_create_dmat(m, n, &A, mem_A.data());
        blasfeo_create_dmat(n, p, &B, mem_B.data());
        blasfeo_create_dmat(m, p, &C, mem_C.data());

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
