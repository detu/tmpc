#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/SizeT.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <boost/graph/adjacency_list.hpp>


namespace tmpc
{
    struct OcpVertex
    {
        OcpVertex() = default;


        OcpVertex(OcpSize const& sz)
        :   size(sz)
        {            
        }


        OcpSize size;
    };


    using OcpSizeGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::bidirectionalS, 
        OcpVertex,
        boost::property<boost::edge_index_t, size_t>
    >;


	inline auto size(OcpSizeGraph const& g)
    {
        return get(&OcpVertex::size, g);
    }
	
	
    /// Property map returning the size of Q matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_Q(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nx(), sz.nx()); }, size_map);
    }


    /// Property map returning the size of R matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_R(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nu(), sz.nu()); }, size_map);
    }


    /// Property map returning the size of S matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_S(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nx(), sz.nu()); }, size_map);
    }


    /// Property map returning the size of C matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_C(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nc(), sz.nx()); }, size_map);
    }


    /// Property map returning the size of D matrix for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_D(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return std::pair(sz.nc(), sz.nu()); }, size_map);
    }


    /// Property map returning the size of ld, ud vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_d(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nc(); }, size_map);
    }
	
	
	/// Property map returning the size of x, lx, ux, q vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_x(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nx(); }, size_map);
    }


    /// Property map returning the size of u, lu, uu, r vectors for a given vertex.
    template <typename SizePropertyMap>
    inline auto size_u(SizePropertyMap size_map)
    {
        return transform_value_property_map(
            [] (OcpSize const& sz) { return sz.nu(); }, size_map);
    }


    /// Property map returning the size of A matrix for a given edge.
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


    /// Property map returning the size of B matrix for a given edge.
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


    /// Property map returning the size of b vector for a given edge.
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
 

    /// Create a tree-structured OcpSizeGraph from two lists:
    /// * the first list defined by the out_gedree iterator is the number out-edges of each node, in breadth-first order.
    /// * the second list defined by the sz iterator is the OcpSize of each node, in breadth-first order.
    ///
    template <typename InIterOutDegree, typename InIterOcpSize>
    inline OcpSizeGraph ocpSizeGraphFromOutDegreeList(InIterOutDegree out_degree, InIterOcpSize sz)
    {
        // Traverse the out-degree list and calculate the total number of vertices.
        auto od = out_degree;
        size_t u = 0;
        size_t v = 1;

        while (u < v)
        {
            v += *od;
            ++u;
            ++od;
        }

        size_t const num_v = u;

        // Create the graph.
        OcpSizeGraph g(num_v);

        // Traverse the out-degree list again, create edges and set node sizes.
        od = out_degree;
        u = 0;
        v = 1;

        while (u < v)
        {
            g[u].size = *sz;

            for (size_t k = 0; k < *od; ++k)
            {
                add_edge(u, v, v - 1 /* edge index */, g);
                ++v;
            }

            ++u;
            ++od;
            ++sz;
        }

        return g;
    }


    /// Create a linear OcpSizeGraph with sized defined by an iterator range [first, last) of OcpSize.
    ///
    template <typename InIterOcpSize>
    inline OcpSizeGraph ocpSizeGraphLinear(InIterOcpSize first, InIterOcpSize last)
    {
        // Create the graph.
        OcpSizeGraph g(std::distance(first, last));

        size_t v = 0;
        for (; first != last; ++first, ++v)
        {
            g[v].size = *first;

            if (v > 0)
                add_edge(v - 1, v, v - 1 /* edge index */, g);
        }

        return g;
    }


    /**
	 * \brief OcpSizeGraph corresponding to a nominal MPC problem with given sizes.
     * 
     * \param nt is the number of control intervals. The total number of nodes is nt + 1.
	 */
	OcpSizeGraph ocpSizeGraphNominalMpc(size_t nt, size_t nx, size_t nu, size_t nc, size_t nct);
}
