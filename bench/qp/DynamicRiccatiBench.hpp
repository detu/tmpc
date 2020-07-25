#pragma once

#include <bench/Benchmark.hpp>


namespace tmpc :: benchmark
{
    inline void riccatiBenchArguments(internal::Benchmark* b) 
    {
        // for (int nx = 1; nx <= 39; ++nx)
        //    b->Args({100, nx, 1});

        for (int nx = 1; nx <= 39; ++nx)
            b->Args({100, nx, nx});
    }
}