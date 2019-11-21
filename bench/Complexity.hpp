#pragma once

#include <cstdlib>
#include <string>
#include <map>
#include <stdexcept>

#include <bench/Benchmark.hpp>


namespace tmpc :: benchmark
{
    using Complexity = std::map<std::string, std::size_t>;

    /// @brief Algorithmic complexity of classical Riccati algorithm
    inline Complexity complexityClassicalRiccati(std::size_t nx, std::size_t nu, std::size_t n)
    {
        using std::pow;

        return {
            // From Frison2013
            {"flops", n * (4 * pow(nx, 3) + 6 * pow(nx, 2) * nu + 3 * nx * pow(nu, 2) + pow(nu, 3) / 3)}
        };
    }


    /// @brief Algorithmic complexity of factorized Riccati algorithm (Frison2013)
    inline Complexity complexityFactorizedRiccati(std::size_t nx, std::size_t nu, std::size_t n)
    {
        using std::pow;

        return {
            // From Frison2013
            {"flops", n * (7 * pow(nx, 3) / 3 + 4 * pow(nx, 2) * nu + 2 * nx * pow(nu, 2) + pow(nu, 3) / 3)}
        };
    }


    template <typename Map>
    inline void setCounters(Map& counters, Complexity const& c)
    {
        for (auto const& v : c)
            counters["f." + v.first] = Counter(v.second, Counter::kIsIterationInvariantRate);

        if (c.find("flops") == c.end())
        {
            std::size_t flops = 0;
            for (auto op : {"add", "mul"})
            {
                auto const it = c.find(op);
                if (it != c.end())
                    flops += it->second;
            }

            counters["flops"] = Counter(flops, Counter::kIsIterationInvariantRate);;
        }
    }
}