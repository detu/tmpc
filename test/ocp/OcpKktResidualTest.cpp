#include <tmpc/ocp/OcpKktResidual.hpp>
#include <tmpc/Testing.hpp>

#include <array>


namespace tmpc :: testing
{
	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	TEST(OcpKktResidualTest, testCtor)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = {2, 1, 1, 0, 0};
		
		std::vector<OcpSize> const sz =
		{
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {1, 0},
			OcpSize {4, 0},
		};

		OcpGraph const g = ocpGraphFromOutDegree(ns.begin());
		auto const size_map = make_iterator_property_map(begin(sz), get(graph::vertex_index, g));
		OcpKktResidual<double> kkt_res {g, size_map};

		for (auto v : graph::vertices(g))
		{
			EXPECT_EQ(get(kkt_res.gx(), v).size(), get(size_map, v).nx());
			EXPECT_EQ(get(kkt_res.gu(), v).size(), get(size_map, v).nu());
		}


		for (auto e : graph::edges(g))
		{
			EXPECT_EQ(get(kkt_res.c(), e).size(), get(size_map, target(e, g)).nx());
		}
	}
}