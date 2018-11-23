#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/ImpactRecorder.hpp>
#include <tmpc/graph/DepthFirstSearch.hpp>

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
		EXPECT_THAT(make_iterator_range(first, last), ElementsAre(
			std::pair(0, 1),
			std::pair(0, 2),
			std::pair(1, 3),
			std::pair(2, 4)
		));
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

		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g)(0), 1);
		EXPECT_EQ(graph::adjacent_vertices(vertex(0, g), g)(1), 2);
		EXPECT_EQ(graph::adjacent_vertices(vertex(1, g), g)(0), 3);
		EXPECT_EQ(graph::adjacent_vertices(vertex(2, g), g)(0), 4);

		// NOTE: use () instead of [] for indexing, because of this:
		/// https://github.com/boostorg/range/issues/83
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(0, g), g)(0)), 0);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(0, g), g)(1)), 1);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(1, g), g)(0)), 2);
		EXPECT_EQ(get(graph::edge_index, g, graph::out_edges(vertex(2, g), g)(0)), 3);
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
			EXPECT_EQ(get(graph::edge_index, g, e), get(graph::vertex_index, g, target(e, g)) - 1);
			++expected_edge;
		}
	}


	TEST(OcpGraphTest, test_parent)
	{
		OcpGraph const g = ocpGraphRobustMpc(3, 2, 1);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_FALSE(parent(0, g).has_value());
		EXPECT_EQ(parent(1, g), 0);
		EXPECT_EQ(parent(2, g), 0);
		EXPECT_EQ(parent(3, g), 1);
		EXPECT_EQ(parent(4, g), 2);
	}


	TEST(OcpGraphTest, test_siblings)
	{
		OcpGraph const g = ocpGraphRobustMpc(3, 2, 1);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_THAT(siblings(0, g), ElementsAre());
		EXPECT_THAT(siblings(1, g), ElementsAre(1, 2));
		EXPECT_THAT(siblings(2, g), ElementsAre(1, 2));
		EXPECT_THAT(siblings(3, g), ElementsAre(3));
		EXPECT_THAT(siblings(4, g), ElementsAre(4));
	}


	TEST(OcpGraphTest, test_impact)
	{
		OcpGraph const g = ocpGraphRobustMpc(3, 2, 1);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		std::vector<size_t> v_impact(num_vertices(g));
		iterator_property_map impact(v_impact.begin(), get(graph::vertex_index, g));
		std::vector<boost::default_color_type> color(num_vertices(g));

		depth_first_visit(g, root(g), graph::dfs_visitor(graph::ImpactRecorder(impact)), 
			make_iterator_property_map(color.begin(), get(graph::vertex_index, g)));
		// graph::recordImpact(g, impact, make_iterator_property_map(color.begin(), get(graph::vertex_index, g)));

		EXPECT_EQ(get(g.impact(), 0), 2);
		EXPECT_EQ(get(g.impact(), 1), 1);
		EXPECT_EQ(get(g.impact(), 2), 1);
		EXPECT_EQ(get(g.impact(), 3), 1);
		EXPECT_EQ(get(g.impact(), 4), 1);

		EXPECT_EQ(get(impact, 0), 2);
		EXPECT_EQ(get(impact, 1), 1);
		EXPECT_EQ(get(impact, 2), 1);
		EXPECT_EQ(get(impact, 3), 1);
		EXPECT_EQ(get(impact, 4), 1);
	}
}