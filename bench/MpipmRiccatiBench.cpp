#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include "RiccatiBench.hpp"


namespace tmpc :: benchmark
{
    // static void BenchmarkArguments(benchmark::internal::Benchmark* b) 
    // {
    //     for (int nx : {2, 5, 10})
    //         for (int nu : {1, 2})
    //             b->Args({100, nx, nu});
    // }


    
    // BENCHMARK(BM_HpipmRiccati)->Apply(BenchmarkArguments);


    BENCHMARK_TEMPLATE(BM_Riccati, MpipmWorkspace<double>)->Args({100, 10, 5});
}
