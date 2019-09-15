#pragma once

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    inline void riccatiBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        for (int nx = 2; nx <= 30; ++nx)
           b->Args({100, nx, 1});

        for (int nx = 2; nx <= 30; ++nx)
            b->Args({100, nx, nx});
    }
}