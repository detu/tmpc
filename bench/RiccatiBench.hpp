#pragma once

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    inline void riccatiBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        for (int nx = 2; nx < 30; ++nx)
            for (int nu = 1; nu < 10; ++nu)
                b->Args({100, nx, nu});
    }
}