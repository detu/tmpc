#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/SizeT.hpp>
#include <tmpc/property_map/PropertyMap.hpp>


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
                if (v > num_intervals)
                    throw std::out_of_range("Vertex index out of range for nominal MPC! (ocpSizeNominalMpc())");

                size_t const actual_nx = v == 0 && first_state_empty ? 0 : nx;
                size_t const actual_nu = v == num_intervals ? 0 : nu;
                size_t const actual_nc = v == num_intervals ? nct : nc;
                
                return OcpSize(actual_nx, actual_nu, actual_nc, ns);                
            }
        );
	}


    /// @brief Property map returning OCP size of robust MPC graph nodes.
	inline auto ocpSizeRobustMpc(OcpGraph const& g, size_t nx, size_t nu, size_t nc = 0, size_t ns = 0, size_t nct = 0, 
        bool first_state_empty = true)
	{
		return make_function_property_map<OcpVertexDescriptor>(
            [&g, nx, nu, nc, ns, nct, first_state_empty] (OcpVertexDescriptor v)
            {
                size_t const actual_nx = in_degree(v, g) == 0 && first_state_empty ? 0 : nx;
                size_t const actual_nu = out_degree(v, g) == 0 ? 0 : nu;
                size_t const actual_nc = out_degree(v, g) == 0 ? nct : nc;

                return OcpSize(actual_nx, actual_nu, actual_nc, ns);
            }
        );
	}


    /// @brief Property map returning OCP size of nominal MHE graph nodes.
	inline auto ocpSizeNominalMhe(size_t num_intervals, size_t nx, size_t nw, size_t nc = 0, size_t ns = 0)
	{
		return make_function_property_map<OcpVertexDescriptor>(
            [num_intervals, nx, nw, nc, ns] (OcpVertexDescriptor v)
            {
                if (v > num_intervals)
                    throw std::out_of_range("Vertex index out of range for nominal MHE! (ocpSizeNominalMhe())");

                return OcpSize(nx, v == num_intervals ? 0 : nw, nc, ns);                
            }
        );
	}
}
