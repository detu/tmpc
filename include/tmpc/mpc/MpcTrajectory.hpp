#pragma once

#include "MpcTrajectoryPoint.hpp"

#include <iterator>

namespace tmpc
{
    template <typename Real>
    using MpcTrajectory = std::vector<MpcTrajectoryPoint<Real>>;

    template <typename Real>
    inline MpcTrajectory<Real> constantMpcTrajectory(MpcTrajectoryPoint<Real> const& p, size_t n)
    {
        MpcTrajectory<Real> tr;
        tr.reserve(n + 1);

        std::fill_n(std::back_inserter(tr), n, p);
        tr.push_back(MpcTrajectoryPoint<Real>(p.x(), DynamicVector<Real>(0)));

        return tr;
    }
}