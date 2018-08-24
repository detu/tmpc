#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    template <typename Kernel, size_t M, size_t N, size_t P>
    static void BM_MTimes(::benchmark::State& state)
    {
        StaticMatrix<Kernel, M, N> A;
        randomize(A);

        StaticMatrix<Kernel, N, P> B;
        randomize(B);

        StaticMatrix<Kernel, M, P> C;
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(C = A * B);
    }


    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 2, 2, 2);
    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 3, 3, 3);
    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 5, 5, 5);
    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 10, 10, 10);
    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 20, 20, 20);
    BENCHMARK_TEMPLATE(BM_MTimes, EigenKernel<double>, 30, 30, 30);
    
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 2, 2, 2);
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 3, 3, 3);
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 5, 5, 5);
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 10, 10, 10);
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 20, 20, 20);
    BENCHMARK_TEMPLATE(BM_MTimes, BlazeKernel<double>, 30, 30, 30);
}
