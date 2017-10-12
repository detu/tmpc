#pragma once

#include <tmpc/ocp/OcpPoint.hpp>

#include <iterator>

namespace tmpc
{
    template <typename Kernel>
    inline auto constantMpcTrajectory(OcpPoint<Kernel> const& p, size_t n)
    {
        std::vector<OcpPoint<Kernel>> tr;
        tr.reserve(n + 1);

        std::fill_n(std::back_inserter(tr), n, p);
        tr.emplace_back(p.x(), StaticVector<Kernel, 0> {});

        return tr;
    }
}