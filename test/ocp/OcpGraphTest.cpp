#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <tmpc/test_tools.hpp>

#include <array>


namespace tmpc :: testing
{
	TEST(OcpGraphTest, test_OutDegreeListToEdgesIterator)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = { 2, 1, 1, 0, 0, };
		using arrayIterator = std::array<size_t, N>::const_iterator;

		tmpc::detail::OutDegreeListToEdgesIterator<arrayIterator> first(ns.begin()), last;

		//for (; first != last; ++first)
		//	std::cout << first->first << ", " << first->second << std::endl;
		EXPECT_EQ(std::distance(first, last), N - 1);
	}


	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	TEST(OcpGraphTest, test_ocpGraphFromOutDegree)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = { 2, 1, 1, 0, 0, };

		OcpGraph const g = ocpGraphFromOutDegree(ns.begin());
		
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		EXPECT_EQ(out_degree(vertex(0, g), g), 2);
		EXPECT_EQ(out_degree(vertex(1, g), g), 1);
		EXPECT_EQ(out_degree(vertex(2, g), g), 1);
		EXPECT_EQ(out_degree(vertex(3, g), g), 0);
		EXPECT_EQ(out_degree(vertex(4, g), g), 0);

		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g)[0], 1);
		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g)[1], 2);
		EXPECT_EQ(graph::adjacent_vertices(vertex(1, g), g)[0], 3);
		EXPECT_EQ(graph::adjacent_vertices(vertex(2, g), g)[0], 4);

		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(0, g), g)[0]), 0);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(0, g), g)[1]), 1);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(1, g), g)[0]), 2);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(2, g), g)[0]), 3);
	}


	// Build a linear graph
    //
    // 0 - 1 - 2 - 3
	TEST(OcpGraphTest, test_ocpGraphLinear)
	{
		size_t const N = 4;
		OcpGraph const g = ocpGraphLinear(N);
	
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		EXPECT_EQ(out_degree(vertex(0, g), g), 1);
		EXPECT_EQ(out_degree(vertex(1, g), g), 1);
		EXPECT_EQ(out_degree(vertex(2, g), g), 1);
		EXPECT_EQ(out_degree(vertex(3, g), g), 0);

		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g).front(), vertex(1, g));
		EXPECT_EQ(graph::adjacent_vertices(vertex(1, g), g).front(), vertex(2, g));
		EXPECT_EQ(graph::adjacent_vertices(vertex(2, g), g).front(), vertex(3, g));

		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g).begin() + 1, graph::adjacent_vertices(vertex(0, g), g).end());
		EXPECT_EQ(graph::adjacent_vertices(vertex(1, g), g).begin() + 1, graph::adjacent_vertices(vertex(1, g), g).end());
		EXPECT_EQ(graph::adjacent_vertices(vertex(2, g), g).begin() + 1, graph::adjacent_vertices(vertex(2, g), g).end());
		EXPECT_EQ(graph::adjacent_vertices(vertex(3, g), g).begin(), graph::adjacent_vertices(vertex(3, g), g).end());

		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(0, g).front()), 0);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(1, g).front()), 1);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(2, g).front()), 2);
	}
}