#include <blasfeo_d_aux.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_s_blas.h>

#include <cblas.h>

#include <benchmark/benchmark.h>

#include <random>
#include <memory>
#include <map>
#include <iostream>


namespace tmpc :: benchmark
{
    using std::size_t;


    static void randomize(size_t m, size_t n, blasfeo_dmat * A)
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(-1.0, 1.0);	

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                BLASFEO_DMATEL(A, i, j) = dis(gen);
    }


    static void randomize(size_t m, size_t n, double * A)
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(-1.0, 1.0);	

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                A[i + j * m] = dis(gen);
    }


    static auto alignedAlloc(size_t bytes, size_t alignment)
    {
        return std::unique_ptr<char[], decltype(&std::free)>(
            reinterpret_cast<char *>(std::aligned_alloc(alignment, bytes)), &std::free);
    }


    // Allocates, holds and and frees memory chunks required by the benchmark.
    struct Mem
    {
        // Disable copying
        Mem(Mem const&) = delete;


        // Move constructor
        Mem(Mem&& rhs)
        :   A_(rhs.A_)
        ,   B_(rhs.B_)
        ,   C_(rhs.C_)
        ,   D_(rhs.D_)
        {
            rhs.A_ = rhs.B_ = rhs.C_ = rhs.D_ = nullptr;
        }


        Mem(size_t m, size_t n, size_t k, size_t alignment)
        :   A_(std::aligned_alloc(alignment, blasfeo_memsize_dmat(k, m)))
        ,   B_(std::aligned_alloc(alignment, blasfeo_memsize_dmat(k, n)))
        ,   C_(std::aligned_alloc(alignment, blasfeo_memsize_dmat(m, n)))
        ,   D_(std::aligned_alloc(alignment, blasfeo_memsize_dmat(m, n)))
        {
        }


        ~Mem()
        {
            free(D_);
            free(C_);
            free(B_);
            free(A_);
        }


        void * A_ = nullptr;
        void * B_ = nullptr;
        void * C_ = nullptr;
        void * D_ = nullptr;
    };


    inline std::ostream& operator<<(std::ostream& os, Mem const& mem)
    {
        return os << "A: " << mem.A_ << "\tB: " << mem.B_ << "\tC: " << mem.C_ << "\tD: " << mem.D_;
    }


    static void BM_gemm_blasfeo(::benchmark::State& state)
    {
        size_t const m = state.range(0);
        size_t const n = state.range(1);
        size_t const k = state.range(2);
        size_t const alignment = state.range(3);

        Mem mem(m, n, k, alignment);

        // std::cout << "Benchmark size " << m << ", " << n << ", " << k << std::endl;
        // std::cout << "Allocated memory pointers " << mem << std::endl;

        blasfeo_dmat A, B, C, D;
        blasfeo_create_dmat(k, m, &A, mem.A_);
        blasfeo_create_dmat(k, n, &B, mem.B_);
        blasfeo_create_dmat(m, n, &C, mem.C_);
        blasfeo_create_dmat(m, n, &D, mem.D_);

        randomize(k, m, &A);
        randomize(k, n, &B);
        randomize(m, n, &C);

        for (auto _ : state)
            blasfeo_dgemm_tn(m, n, k, 
                1.0, 
                &A, 0, 0, 
                &B, 0, 0, 
                1.0, 
                &C, 0, 0, 
                &D, 0, 0);
    }


    static void BM_gemm_blasfeo_reuse_memory(::benchmark::State& state)
    {
        size_t const m = state.range(0);
        size_t const n = state.range(1);
        size_t const k = state.range(2);
        size_t const alignment = state.range(3);

        // Data type describing benchmark settings
        using Settings = std::array<size_t, 4>;
        Settings const settings {m, n, k, alignment};

        // Map from settings to memory chunks
        // (persistent between function calls due to static)
        static std::map<Settings, Mem> mem_map;

        // Check if we have already allocated memory for these settings
        auto mem = mem_map.find(settings);
        if (mem == mem_map.end())
            // Allocate memory chunks if not already allocated
            mem = mem_map.emplace(settings, Mem(m, n, k, alignment)).first;

        // std::cout << "Benchmark size " << m << ", " << n << ", " << k << std::endl;
        // std::cout << "Allocated memory pointers " << mem << std::endl;

        blasfeo_dmat A, B, C, D;
        blasfeo_create_dmat(k, m, &A, mem->second.A_);
        blasfeo_create_dmat(k, n, &B, mem->second.B_);
        blasfeo_create_dmat(m, n, &C, mem->second.C_);
        blasfeo_create_dmat(m, n, &D, mem->second.D_);

        randomize(k, m, &A);
        randomize(k, n, &B);
        randomize(m, n, &C);

        for (auto _ : state)
            blasfeo_dgemm_tn(m, n, k, 
                1.0, 
                &A, 0, 0, 
                &B, 0, 0, 
                1.0, 
                &C, 0, 0, 
                &D, 0, 0);
    }


    static void BM_gemm_cblas(::benchmark::State& state)
    {
        size_t const m = state.range(0);
        size_t const n = state.range(1);
        size_t const k = state.range(2);

        auto A = std::make_unique<double []>(k * m);
        auto B = std::make_unique<double []>(k * n);
        auto C = std::make_unique<double []>(m * n);

        randomize(k, m, A.get());
        randomize(k, n, B.get());
        randomize(m, n, C.get());

        for (auto _ : state)
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k,
		    1.0, A.get(), k, B.get(), k, 1.0, C.get(), m);

    }

    
    BENCHMARK(BM_gemm_blasfeo)
        ->Args({2, 2, 2, 0x20})
        ->Args({3, 3, 3, 0x20})
        ->Args({5, 5, 5, 0x20})
        ->Args({10, 10, 10, 0x20})
        ->Args({20, 20, 20, 0x20})
        ->Args({30, 30, 30, 0x20})
        ->Args({2, 2, 2, 0x1000})
        ->Args({3, 3, 3, 0x1000})
        ->Args({5, 5, 5, 0x1000})
        ->Args({10, 10, 10, 0x1000})
        ->Args({20, 20, 20, 0x1000})
        ->Args({30, 30, 30, 0x1000});


#if 0
    // Run benchmarks in normal order
    BENCHMARK(BM_gemm_blasfeo_reuse_memory)
        ->Args({2, 2, 2, 0x20})
        ->Args({3, 3, 3, 0x20})
        ->Args({5, 5, 5, 0x20})
        ->Args({10, 10, 10, 0x20})
        ->Args({20, 20, 20, 0x20})
        ->Args({30, 30, 30, 0x20});
#else
    // Run benchmarks in reverse order
    BENCHMARK(BM_gemm_blasfeo_reuse_memory)
        ->Args({30, 30, 30, 0x20})
        ->Args({20, 20, 20, 0x20})
        ->Args({10, 10, 10, 0x20})
        ->Args({5, 5, 5, 0x20})
        ->Args({3, 3, 3, 0x20})
        ->Args({2, 2, 2, 0x20});
#endif


    BENCHMARK(BM_gemm_cblas)
        ->Args({2, 2, 2})
        ->Args({3, 3, 3})
        ->Args({5, 5, 5})
        ->Args({10, 10, 10})
        ->Args({20, 20, 20})
        ->Args({30, 30, 30});
}
