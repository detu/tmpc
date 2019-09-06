#include <tmpc/math/Llh.hpp>

#include <blaze/Math.h>
#include <Eigen/Dense>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real>
    static void BM_Cholesky_Blaze_potrf(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(N, N);        
        makePositiveDefinite(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);        
        for (auto _ : state)
        {
            B = A;
            potrf(B, 'L');
        }
    }


    template <typename Real>
    static void BM_Cholesky_Blaze_llh_Dynamic(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(N, N);        
        makePositiveDefinite(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);        
        for (auto _ : state)
            blaze::llh(A, B);
    }


    template <typename Real>
    static void BM_Cholesky_Blaze_llh_SymmetricDynamic(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> A(N);        
        makePositiveDefinite(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);        
        for (auto _ : state)
            blaze::llh(A, B);
    }


    template <typename Real>
    static void BM_Cholesky_Blaze_llh_SymmetricToLowerDynamic(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> A(N);        
        makePositiveDefinite(A);

        blaze::LowerMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> B(N);        
        for (auto _ : state)
            blaze::llh(A, B);
    }


    template <typename Real>
    static void BM_Cholesky_Blaze_llh_DeclsymToDynamic(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(N, N);        
        makePositiveDefinite(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);        
        for (auto _ : state)
            blaze::llh(declsym(A), B);
    }


    template <typename Real>
    static void BM_Cholesky_Eigen_LLT_Dynamic(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
        size_t const N = state.range(0);
        
        // Make a positive definite matrix A
        Matrix A = Matrix::Random(N, N);
        A = A.transpose() * A;

        Eigen::LLT<Matrix> llt(N);
        for (auto _ : state)
            ::benchmark::DoNotOptimize(llt.compute(A));
    }


    template <typename Real>
    static void BM_Cholesky_Eigen_LLT_SelfadjointDynamic(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
        size_t const N = state.range(0);
        
        // Make a positive definite matrix A
        Matrix A = Matrix::Random(N, N);
        A = A.transpose() * A;

        Eigen::LLT<Matrix> llt(N);
        for (auto _ : state)
            ::benchmark::DoNotOptimize(llt.compute(A.template selfadjointView<Eigen::Lower>()));
    }


    template <typename Real, int N>
    static void BM_Cholesky_Eigen_LLT_Static(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, N, N>;
        
        // Make a positive definite matrix A
        Matrix A = Matrix::Random();
        A = A.transpose() * A;

        Eigen::LLT<Matrix> llt(N);
        for (auto _ : state)
            ::benchmark::DoNotOptimize(llt.compute(A));
    }


    template <typename Real, int N, bool SO1, bool SO2>
    static void BM_Cholesky_tmpc_llh_Static(::benchmark::State& state)
    {
        // Make a positive definite matrix A
        blaze::SymmetricMatrix<blaze::StaticMatrix<Real, N, N, SO1>> A;
        makePositiveDefinite(A);

        blaze::LowerMatrix<blaze::StaticMatrix<Real, N, N, SO2>> L;

        for (auto _ : state)
            tmpc::llh(A, L);
    }


    static void choleskyBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        b->Arg(1)->Arg(2)->Arg(5)->Arg(10)->Arg(35);
    }

    
    BENCHMARK_TEMPLATE(BM_Cholesky_Blaze_llh_Dynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Blaze_llh_SymmetricDynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Blaze_llh_SymmetricToLowerDynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Blaze_llh_DeclsymToDynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Blaze_potrf, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Dynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_SelfadjointDynamic, double)->Apply(choleskyBenchArguments);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Static, double, 1);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Static, double, 2);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Static, double, 5);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Static, double, 10);
    BENCHMARK_TEMPLATE(BM_Cholesky_Eigen_LLT_Static, double, 35);
    
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::rowMajor, blaze::rowMajor);
    
}
