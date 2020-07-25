#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/Testing.hpp>

#include <array>


namespace tmpc :: testing
{
	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	TEST(OcpTreeTest, testFromOutDegree)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = { 2, 1, 1, 0, 0, };

		OcpTree const g {ns};
		
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		EXPECT_EQ(out_degree(vertex(0, g), g), 2);
		EXPECT_EQ(out_degree(vertex(1, g), g), 1);
		EXPECT_EQ(out_degree(vertex(2, g), g), 1);
		EXPECT_EQ(out_degree(vertex(3, g), g), 0);
		EXPECT_EQ(out_degree(vertex(4, g), g), 0);

		EXPECT_EQ(children(vertex(0, g), g)[0], 1);
		EXPECT_EQ(children(vertex(0, g), g)[1], 2);
		EXPECT_EQ(children(vertex(1, g), g)[0], 3);
		EXPECT_EQ(children(vertex(2, g), g)[0], 4);

		// NOTE: use () instead of [] for indexing, because of this:
		/// https://github.com/boostorg/range/issues/83
		EXPECT_EQ(out_edges(vertex(0, g), g)[0], 0);
		EXPECT_EQ(out_edges(vertex(0, g), g)[1], 1);
		EXPECT_EQ(out_edges(vertex(1, g), g)[0], 2);
		EXPECT_EQ(out_edges(vertex(2, g), g)[0], 3);
	}


	// Build a linear graph
    //
    // 0 - 1 - 2 - 3
	TEST(OcpTreeTest, testLinear)
	{
		size_t const N = 4;
		OcpTree const g(N);
	
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		for (auto v : vertices(g))
		{
			auto const nc_expected = v + 1 < N ? 1 : 0;
			EXPECT_EQ(out_degree(v, g), nc_expected);
			EXPECT_EQ(std::size(children(v, g)), nc_expected);

			if (nc_expected > 0)
				EXPECT_EQ(children(v, g).front(), v + 1);
			else
				EXPECT_TRUE(std::empty(children(v, g)));
		}
	}


	// Build a robust MPC tree
	TEST(OcpTreeTest, testRobustMpc)
	{
		OcpTree const g(4, 2, 2);

		std::vector<std::pair<size_t, size_t>> expected_edges = {
			{0, 1}, {0, 2},
			{1, 3}, {1, 4}, {2, 5}, {2, 6},
			{3, 7}, {4, 8}, {5, 9}, {6, 10}
		};
	
		ASSERT_EQ(num_vertices(g), expected_edges.size() + 1);
		ASSERT_EQ(num_edges(g), expected_edges.size());

		auto expected_edge = expected_edges.begin();
		for (auto e : edges(g))
		{
			EXPECT_EQ(source(e, g), expected_edge->first);
			EXPECT_EQ(target(e, g), expected_edge->second);
			EXPECT_EQ(e, target(e, g) - 1);
			++expected_edge;
		}
	}


	// Check that the vertices(g) returns a correct sequence of vertices.
	//
	TEST(OcpTreeTest, testVertices)
	{
		size_t const N = 4;
		OcpTree const g(N);
		auto const v = vertices(g);

		ASSERT_EQ(std::size(v), 4);
		EXPECT_EQ(v[0], 0);
		EXPECT_EQ(v[1], 1);
		EXPECT_EQ(v[2], 2);
		EXPECT_EQ(v[3], 3);
	}


	// Check that the vertices(g) | std::views::reverse returns a correct sequence of vertices.
	//
	TEST(OcpTreeTest, testReverseVertices)
	{
		size_t const N = 4;
		OcpTree const g(N);
		auto const v = vertices(g) | std::views::reverse;

		ASSERT_EQ(std::size(v), 4);
		EXPECT_EQ(v[0], 3);
		EXPECT_EQ(v[1], 2);
		EXPECT_EQ(v[2], 1);
		EXPECT_EQ(v[3], 0);
	}
	
	
	TEST(OcpTreeTest, testParent)
	{
		OcpTree const g(3, 1, 2);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_FALSE(parent(0, g).has_value());
		EXPECT_EQ(parent(1, g), 0);
		EXPECT_EQ(parent(2, g), 0);
		EXPECT_EQ(parent(3, g), 1);
		EXPECT_EQ(parent(4, g), 2);
	}


	TEST(OcpTreeTest, testSiblings)
	{
		OcpTree const g(3, 1, 2);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_EQ(std::size(siblings(0, g)), 1);
		EXPECT_EQ(siblings(0, g)[0], 0);
		
		EXPECT_EQ(std::size(siblings(1, g)), 2);
		EXPECT_EQ(siblings(1, g)[0], 1);
		EXPECT_EQ(siblings(1, g)[1], 2);
		
		EXPECT_EQ(std::size(siblings(2, g)), 2);
		EXPECT_THAT(siblings(2, g)[0], 1);
		EXPECT_THAT(siblings(2, g)[1], 2);

		EXPECT_EQ(std::size(siblings(3, g)), 1);
		EXPECT_EQ(siblings(3, g)[0], 3);

		EXPECT_EQ(std::size(siblings(4, g)), 1);
		EXPECT_EQ(siblings(4, g)[0], 4);
	}


	TEST(OcpTreeTest, testScenarioCount)
	{
		OcpTree const g(3, 1, 2);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_EQ(g.scenarioCount(0), 2);
		EXPECT_EQ(g.scenarioCount(1), 1);
		EXPECT_EQ(g.scenarioCount(2), 1);
		EXPECT_EQ(g.scenarioCount(3), 1);
		EXPECT_EQ(g.scenarioCount(4), 1);
	}


	TEST(OcpTreeTest, testDepth)
	{
		OcpTree const g(3, 1, 2);

		ASSERT_EQ(num_vertices(g), 5);
		ASSERT_EQ(num_edges(g), 4);

		EXPECT_EQ(g.depth(0), 0);
		EXPECT_EQ(g.depth(1), 1);
		EXPECT_EQ(g.depth(2), 1);
		EXPECT_EQ(g.depth(3), 2);
		EXPECT_EQ(g.depth(4), 2);
	}


	TEST(OcpTreeTest, testEqual)
	{
		OcpTree const g1 {2, 1, 1, 0, 0};
		OcpTree const g2 {2, 1, 1, 0, 0};
		OcpTree const g3 {2, 1, 2, 0, 0, 0};

		EXPECT_EQ(g1, g1);
		EXPECT_EQ(g1, g2);
		EXPECT_EQ(g2, g1);
		EXPECT_NE(g1, g3);
		EXPECT_NE(g3, g1);

		EXPECT_EQ(g2, g2);
		EXPECT_NE(g2, g3);
		EXPECT_NE(g3, g2);

		EXPECT_EQ(g3, g3);
	}


	TEST(OcpTreeTest, testCopyCtor)
	{
		OcpTree const g1 {2, 1, 1, 0, 0};
		OcpTree const g2 {g1};

		EXPECT_EQ(g2, g1);
	}


	TEST(OcpTreeTest, testMoveCtor)
	{
		OcpTree const g1 {2, 1, 1, 0, 0};
		OcpTree g2 {g1};
		OcpTree g3 {std::move(g2)};

		EXPECT_EQ(g3, g1);
	}


	TEST(OcpTreeTest, testLeafVertices)
	{
		OcpTree const g {2, 1, 1, 0, 0};
		auto const leaves = g.leafVertices();

		EXPECT_EQ(leaves[0], 3);
		EXPECT_EQ(leaves[1], 4);
	}


	TEST(OcpTreeTest, testBranchVertices)
	{
		OcpTree const g {2, 1, 1, 0, 0};
		auto const branches = g.branchVertices();

		EXPECT_EQ(branches[0], 0);
		EXPECT_EQ(branches[1], 1);
		EXPECT_EQ(branches[2], 2);
	}
}