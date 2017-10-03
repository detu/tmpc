#pragma once

#include <tmpc/ocp/OcpPoint.hpp>

#include <ostream>


namespace tmpc
{
    template <typename Kernel>
    inline std::ostream& operator<<(std::ostream& os, OcpPoint<Kernel> const& p)
    {
        os << "x=" << trans(p.x()) << "\tu=" << trans(p.u());
    }
}

