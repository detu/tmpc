#pragma once

#include <tmpc/ocp/OcpTree.hpp>


namespace tmpc
{
    /// @brief Defines the tree-structured OCP size.
    ///
    template <typename Size>
    concept OcpSize = requires(Size size, OcpVertex v)
    {
        size.graph();
        size.nx(v);
        size.nu(v);
        size.nc(v);
    };


    template <OcpSize Size1, OcpSize Size2>
    inline bool operator==(Size1 const& a, Size2 const& b)
	{
		if (a.graph() != b.graph())
			return false;

		for (auto v : vertices(a.graph()))
			if (a.nx(v) != b.nx(v) || a.nu(v) != b.nu(v) || a.nc(v) != b.nc(v) || a.ns(v) != b.ns(v))
				return false;

		return true;
	}

	
    template <OcpSize Size1, OcpSize Size2>
	inline bool operator!=(Size1 const& a, Size2 const& b)
	{
		return !(a == b);
	}


    /// @brief Number of states at the source vertex of e.
    template <OcpSize Size>
    inline auto sourceNx(Size const& size, OcpEdge e)
    {
        return size.nx(source(e, size.graph()));
    }


    /// @brief Number of controls at the source vertex of e.
    template <OcpSize Size>
    inline auto sourceNu(Size const& size, OcpEdge e)
    {
        return size.nu(source(e, size.graph()));
    }


    /// @brief Number of states at the target vertex of e.
    template <OcpSize Size>
    inline auto targetNx(Size const& size, OcpEdge e)
    {
        return size.nx(target(e, size.graph()));
    }
}