#include "RiccatiBench.hpp"


namespace tmpc :: benchmark
{
    void riccatiBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        for (int nx = 2; nx <= 30; ++nx)
            //for (int nu = 1; nu <= 10; ++nu)
            for (int nu : {1, 5})
                b->Args({100, nx, nu});
    }
}