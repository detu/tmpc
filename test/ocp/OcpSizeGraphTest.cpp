#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <tmpc/test_tools.hpp>

#include <array>


namespace tmpc :: testing
{
	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	TEST(OcpSizeGraphTest, test_OcpSizeGraphFromOutDegreeList)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = { 2, 1, 1, 0, 0, };
		std::array<OcpSize, N> const sz = {
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {1, 0},
			OcpSize {4, 0},
		};

		OcpSizeGraph const g = ocpSizeGraphFromOutDegreeList(ns.begin(), sz.begin());
	
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		for (size_t i = 0; i < N; ++i)
			EXPECT_EQ(g[i].size, sz[i]);

		EXPECT_EQ(out_degree(0, g), 2);
		EXPECT_EQ(out_degree(1, g), 1);
		EXPECT_EQ(out_degree(2, g), 1);
		EXPECT_EQ(out_degree(3, g), 0);
		EXPECT_EQ(out_degree(4, g), 0);

		EXPECT_EQ(adjacent_vertices(0, g).first[0], 1);
		EXPECT_EQ(adjacent_vertices(0, g).first[1], 2);
		EXPECT_EQ(adjacent_vertices(1, g).first[0], 3);
		EXPECT_EQ(adjacent_vertices(2, g).first[0], 4);
	}


	// Build a linear graph
    //
    // 0 - 1 - 2 - 3
	TEST(OcpSizeGraphTest, test_OcpSizeGraphLinear)
	{
		size_t const N = 4;
		std::array<OcpSize, N> const sz = {
			OcpSize {2, 1},
			OcpSize {3, 4},
			OcpSize {5, 6},
			OcpSize {7, 8},
		};

		OcpSizeGraph const g = ocpSizeGraphLinear(sz.begin(), sz.end());
	
		EXPECT_EQ(num_vertices(g), N);
		EXPECT_EQ(num_edges(g), N - 1);

		for (size_t i = 0; i < N; ++i)
			EXPECT_EQ(g[i].size, sz[i]);

		EXPECT_EQ(out_degree(0, g), 1);
		EXPECT_EQ(out_degree(1, g), 1);
		EXPECT_EQ(out_degree(2, g), 1);
		EXPECT_EQ(out_degree(3, g), 0);

		EXPECT_EQ(*adjacent_vertices(0, g).first, 1);
		EXPECT_EQ(*adjacent_vertices(1, g).first, 2);
		EXPECT_EQ(*adjacent_vertices(2, g).first, 3);

		EXPECT_EQ(adjacent_vertices(0, g).first + 1, adjacent_vertices(0, g).second);
		EXPECT_EQ(adjacent_vertices(1, g).first + 1, adjacent_vertices(1, g).second);
		EXPECT_EQ(adjacent_vertices(2, g).first + 1, adjacent_vertices(2, g).second);
		EXPECT_EQ(adjacent_vertices(3, g).first, adjacent_vertices(3, g).second);
	}
}