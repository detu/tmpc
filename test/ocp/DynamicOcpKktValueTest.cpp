#include <tmpc/ocp/DynamicOcpKktValue.hpp>
#include <tmpc/Testing.hpp>

#include <array>


namespace tmpc :: testing
{
	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	TEST(DynamicOcpKktValueTest, testCtor)
	{
		size_t const N = 5;
		std::array<size_t, N> const ns = {2, 1, 1, 0, 0};
		
		OcpTree const g {ns};
		DynamicOcpSize const sz {g, {
			{2, 2},
			{2, 2},
			{2, 2},
			{1, 0},
			{4, 0},
		}};
		
		DynamicOcpKktValue<double> kkt_res {sz};

		for (auto v : vertices(g))
		{
			EXPECT_EQ(kkt_res.gx(v).size(), sz.nx(v));
			EXPECT_EQ(kkt_res.gu(v).size(), sz.nu(v));
		}


		for (auto e : edges(g))
		{
			EXPECT_EQ(kkt_res.c(e).size(), sz.nx(target(e, sz.graph())));
		}
	}
}