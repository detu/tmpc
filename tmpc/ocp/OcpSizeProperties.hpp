#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/SizeT.hpp>
#include <tmpc/core/PropertyMap.hpp>


namespace tmpc
{
    /// @brief Property map returning the size of Q matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_Q(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nx(), sz.nx()); }, size_map);
    }


    /// @brief Property map returning the size of R matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_R(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nu(), sz.nu()); }, size_map);
    }


    /// @brief Property map returning the size of S matrix for a given vertex.
    ///
    template <typename SizePropertyMap>
    inline auto size_S(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nu(), sz.nx()); }, size_map);
    }


    /// @brief Property map returning the size of C matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_C(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nc(), sz.nx()); }, size_map);
    }


    /// @brief Property map returning the size of D matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_D(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nc(), sz.nu()); }, size_map);
    }


    /// @brief Property map returning the size of ld, ud vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_d(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nc(); }, size_map);
    }
	
	
	/// @brief Property map returning the size of x, lx, ux, q vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_x(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nx(); }, size_map);
    }


    /// @brief Property map returning the size of u, lu, uu, r vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_u(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nu(); }, size_map);
    }


    /// @brief Property map returning the size of A matrix for a given edge.
    template <typename SizePropertyMap>
    inline auto size_A(SizePropertyMap size_map, OcpGraph const& g)
    {
        return make_function_property_map<OcpEdgeDescriptor>(
            [size_map, &g] (OcpEdgeDescriptor e) 
            { 
                return std::pair(
                    get(size_map, target(e, g)).nx(), 
                    get(size_map, source(e, g)).nx()); 
            }
        );
    }


    /// @brief Property map returning the size of B matrix for a given edge.
    template <typename SizePropertyMap>
    inline auto size_B(SizePropertyMap size_map, OcpGraph const& g)
    {
        return make_function_property_map<OcpEdgeDescriptor>(
            [size_map, &g] (OcpEdgeDescriptor e) 
            { 
                return std::pair(
                    get(size_map, target(e, g)).nx(), 
                    get(size_map, source(e, g)).nu()); 
            }
        );
    }


    /// @brief Property map returning the size of b vector for a given edge.
    template <typename SizePropertyMap>
    inline auto size_b(SizePropertyMap size_map, OcpGraph const& g)
    {
        return make_function_property_map<OcpEdgeDescriptor>(
            [size_map, &g] (OcpEdgeDescriptor e) 
            { 
                return get(size_map, target(e, g)).nx(); 
            }
        );
    }


	/// @brief Property map returning OCP size of nominal MPC graph nodes.
	inline auto ocpSizeNominalMpc(size_t num_intervals, size_t nx, size_t nu, size_t nc = 0, size_t ns = 0, size_t nct = 0, 
        bool first_state_empty = true)
	{
		return make_function_property_map<OcpVertexDescriptor>(
            [num_intervals, nx, nu, nc, ns, nct, first_state_empty] (OcpVertexDescriptor v)
            {
                if (v == 0)
                    return OcpSize(first_state_empty ? 0 : nx, nu, nc, ns);
                if (v < num_intervals)
                    return OcpSize(nx, nu, nc, ns);
                if (v == num_intervals)
                    return OcpSize(nx, 0, nct, ns);

                throw std::out_of_range("Vertex index out of range for nominal MPC! (ocpSizeNominalMpc())");
            }
        );
	}
}
