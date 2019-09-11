#include <blaze/Math.h>
#include <Eigen/Dense>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real>
    static void BM_syrk_Blaze_Dynamic(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        size_t const N = state.range(1);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(M, N);        
        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);

        randomize(A);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = trans(A) * A);
    }


    template <typename Real>
    static void BM_syrk_Blaze_SymmetricDynamic(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        size_t const N = state.range(1);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(M, N);        
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> B(N);

        randomize(A);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = declsym(trans(A) * A));
    }


    template <typename Real, size_t M, size_t N>
    static void BM_syrk_Blaze_Static(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;        
        blaze::StaticMatrix<Real, N, N, blaze::columnMajor> B;

        randomize(A);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = trans(A) * A);
    }


    template <typename Real, size_t M, size_t N>
    static void BM_syrk_Blaze_Static_Declsym(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;        
        blaze::StaticMatrix<Real, N, N, blaze::columnMajor> B;

        randomize(A);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = declsym(trans(A) * A));
    }


    template <typename Real, size_t M, size_t N>
    static void BM_syrk_Blaze_Static_Symmetric_Declsym(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;
        blaze::SymmetricMatrix<blaze::StaticMatrix<Real, N, N, blaze::columnMajor>> B;

        randomize(A);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = declsym(trans(A) * A));
    }


    template <typename Real, size_t M, size_t N, bool SO1, bool SO2>
    void syrkStaticLoop(blaze::StaticMatrix<Real, M, N, SO1> const& A, blaze::StaticMatrix<Real, N, N, SO2>& B)
    {
        for (size_t j = 0; j < N; ++j)
        {
            for (size_t i = j; i < N; ++i)
                B(i, j) = dot(column(A, i), column(A, j));
        }
    }


    template <typename Real, size_t M, size_t N>
    static void BM_syrk_Blaze_Static_Loop(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;        
        blaze::StaticMatrix<Real, N, N, blaze::columnMajor> B;

        randomize(A);
        
        for (auto _ : state)
        {
            syrkStaticLoop(A, B);
        }
    }


    template <typename Real>
    static void BM_syrk_Eigen_Dynamic(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

        size_t const M = state.range(0);
        size_t const N = state.range(1);
        Matrix const A = Matrix::Random(M, N);        
        Matrix B(N, N);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = A.transpose() * A);
    }


    template <typename Real>
    static void BM_syrk_Eigen_RankUpdateDynamic(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

        size_t const M = state.range(0);
        size_t const N = state.range(1);
        Matrix const A = Matrix::Random(M, N);        
        Matrix B(N, N);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B.template selfadjointView<Eigen::Lower>().rankUpdate(A.transpose()));
    }


    template <typename Real>
    static void BM_syrk_Eigen_TriangularDynamic(::benchmark::State& state)
    {
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

        size_t const M = state.range(0);
        size_t const N = state.range(1);
        Matrix const A = Matrix::Random(M, N);        
        Matrix B(N, N);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B.template triangularView<Eigen::Lower>() = A.transpose() * A);
    }


    template <typename Real, int M, int N>
    static void BM_syrk_Eigen_Static(::benchmark::State& state)
    {
        using MatrixA = Eigen::Matrix<Real, M, N>;
        using MatrixB = Eigen::Matrix<Real, N, N>;

        MatrixA const A = MatrixA::Random();
        MatrixB B;
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B = A.transpose() * A);
    }


    template <typename Real, int M, int N>
    static void BM_syrk_Eigen_RankUpdateStatic(::benchmark::State& state)
    {
        using MatrixA = Eigen::Matrix<Real, M, N>;
        using MatrixB = Eigen::Matrix<Real, N, N>;

        MatrixA const A = MatrixA::Random();
        MatrixB B;
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B.template selfadjointView<Eigen::Lower>().rankUpdate(A.transpose()));
    }


    template <typename Real, int M, int N>
    static void BM_syrk_Eigen_TriangularStatic(::benchmark::State& state)
    {
        using MatrixA = Eigen::Matrix<Real, M, N>;
        using MatrixB = Eigen::Matrix<Real, N, N>;

        MatrixA const A = MatrixA::Random();
        MatrixB B;
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(B.template triangularView<Eigen::Lower>() = A.transpose() * A);
    }


#if BLAZE_BLAS_MODE
    template <typename Real>
    static void BM_syrk_Blaze_cblas(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        size_t const N = state.range(1);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(M, N);        
        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);

        randomize(A);
        B = 0.;
        
        for (auto _ : state)
            cblas_dsyrk(CblasColMajor, CblasLower, CblasTrans, 
                N, M, 1., data(A), spacing(A), 0., data(B), spacing(B));
    }
#endif


    static void syrkBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        b->Args({4, 5})->Args({30, 35});
    }


    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Dynamic, double)->Apply(syrkBenchArguments);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_SymmetricDynamic, double)->Apply(syrkBenchArguments);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static, double, 30, 35);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Declsym, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Declsym, double, 30, 35);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Symmetric_Declsym, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Symmetric_Declsym, double, 30, 35);

    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Loop, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_Static_Loop, double, 30, 35);

    BENCHMARK_TEMPLATE(BM_syrk_Eigen_Dynamic, double)->Apply(syrkBenchArguments);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_RankUpdateDynamic, double)->Apply(syrkBenchArguments);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_Static, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_Static, double, 30, 35);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_RankUpdateStatic, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_RankUpdateStatic, double, 30, 35);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_TriangularDynamic, double)->Apply(syrkBenchArguments);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_TriangularStatic, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_syrk_Eigen_TriangularStatic, double, 30, 35);

#if BLAZE_BLAS_MODE
    BENCHMARK_TEMPLATE(BM_syrk_Blaze_cblas, double)->Apply(syrkBenchArguments);
#endif
}
