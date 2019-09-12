#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/Testing.hpp>

#include <array>
#include <iostream>


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

		EXPECT_EQ(get(graph::edge_index, g, *next(graph::out_edges(vertex(0, g), g).begin(), 0)), 0);
		EXPECT_EQ(get(graph::edge_index, g, *next(graph::out_edges(vertex(0, g), g).begin(), 1)), 1);
		EXPECT_EQ(get(graph::edge_index, g, *next(graph::out_edges(vertex(1, g), g).begin(), 0)), 2);
		EXPECT_EQ(get(graph::edge_index, g, *next(graph::out_edges(vertex(2, g), g).begin(), 0)), 3);
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


	// Build a robust MPC tree
	TEST(OcpGraphTest, test_ocpGraphRobustMpc)
	{
		OcpGraph const g = ocpGraphRobustMpc(4, 2, 2);

		std::vector<std::pair<size_t, size_t>> expected_edges = {
			{0, 1}, {0, 2},
			{1, 3}, {1, 4}, {2, 5}, {2, 6},
			{3, 7}, {4, 8}, {5, 9}, {6, 10}
		};
	
		ASSERT_EQ(num_vertices(g), expected_edges.size() + 1);
		ASSERT_EQ(num_edges(g), expected_edges.size());

		auto expected_edge = expected_edges.begin();
		for (auto e : graph::edges(g))
		{
			EXPECT_EQ(source(e, g), expected_edge->first);
			EXPECT_EQ(target(e, g), expected_edge->second);
			++expected_edge;
		}
	}


	// Check that the vertices(g) returns a correct sequence of vertices.
	//
	TEST(OcpGraphTest, test_vertices)
	{
		size_t const N = 4;
		OcpGraph const g = ocpGraphLinear(N);

		EXPECT_THAT(graph::vertices(g), ElementsAreArray({0, 1, 2, 3}));
	}


	// Check that the reverse(vertices(g)) returns a correct sequence of vertices.
	//
	// NOTE: It does not work however because of incorrect behavior of reversed_range:
	// https://github.com/boostorg/range/issues/82
	//
	// The test is disabled until the issue is fixed.
	TEST(OcpGraphTest, DISABLED_test_reverseVertices)
	{
		size_t const N = 4;
		OcpGraph const g = ocpGraphLinear(N);

		EXPECT_THAT(reverse(graph::vertices(g)), ElementsAreArray({3, 2, 1, 0}));
	}
}