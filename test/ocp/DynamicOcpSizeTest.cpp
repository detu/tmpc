#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(DynamicOcpSizeTest, testInitializerListCtor)
	{
		DynamicOcpSize const s {
			{3, 2, 1},
			{5, 4},
			{7, 1},
			{1, 0},
			{4, 1},
		};

		EXPECT_EQ(s.nx(0), 3);
		EXPECT_EQ(s.nu(0), 2);
		EXPECT_EQ(s.nc(0), 1);
		
		EXPECT_EQ(s.nx(1), 5);
		EXPECT_EQ(s.nu(1), 4);
		EXPECT_EQ(s.nc(1), 0);

		EXPECT_EQ(s.nx(2), 7);
		EXPECT_EQ(s.nu(2), 1);
		EXPECT_EQ(s.nc(2), 0);

		EXPECT_EQ(s.nx(3), 1);
		EXPECT_EQ(s.nu(3), 0);
		EXPECT_EQ(s.nc(3), 0);

		EXPECT_EQ(s.nx(4), 4);
		EXPECT_EQ(s.nu(4), 1);
		EXPECT_EQ(s.nc(4), 0);
	}


	TEST(DynamicOcpSizeTest, testGraphAndInitializerListCtor)
	{
		std::vector<size_t> const ns = { 2, 1, 1, 0, 0, };
		
		DynamicOcpSize const s {OcpTree {ns}, {
			{3, 2, 1},
			{5, 4},
			{7, 1},
			{1, 0},
			{4, 1},
		}};

		EXPECT_EQ(s.nx(0), 3);
		EXPECT_EQ(s.nu(0), 2);
		EXPECT_EQ(s.nc(0), 1);
		
		EXPECT_EQ(s.nx(1), 5);
		EXPECT_EQ(s.nu(1), 4);
		EXPECT_EQ(s.nc(1), 0);

		EXPECT_EQ(s.nx(2), 7);
		EXPECT_EQ(s.nu(2), 1);
		EXPECT_EQ(s.nc(2), 0);

		EXPECT_EQ(s.nx(3), 1);
		EXPECT_EQ(s.nu(3), 0);
		EXPECT_EQ(s.nc(3), 0);

		EXPECT_EQ(s.nx(4), 4);
		EXPECT_EQ(s.nu(4), 1);
		EXPECT_EQ(s.nc(4), 0);
	}


	TEST(DynamicOcpSizeTest, testFunctorCtor)
	{
		std::vector<size_t> const ns = { 2, 1, 1, 0, 0, };
		std::vector<OcpVertexSize> const vs {
			{3, 2, 1},
			{5, 4},
			{7, 1},
			{1, 0},
			{4, 1}
		};
		
		DynamicOcpSize const s {OcpTree {ns},
			[&vs] (OcpVertex v, OcpTree const& g)
			{
				return vs[v];
			},
		};

		for (auto v : vertices(s.graph()))
		{
			EXPECT_EQ(s.nx(v), vs[v].nx_);
			EXPECT_EQ(s.nu(v), vs[v].nu_);
			EXPECT_EQ(s.nc(v), vs[v].nc_);
		}
	}


	TEST(DynamicOcpSizeTest, testCopyCtor)
	{
		DynamicOcpSize const s1 {
			{3, 2, 1},
			{5, 4},
			{7, 1},
			{1, 0},
			{4, 1},
		};


		DynamicOcpSize ss1 = s1;
		EXPECT_EQ(ss1, s1);
	}


	TEST(DynamicOcpSizeTest, testNominalMpc)
	{
		OcpTree const g(4);

		size_t const nx = 5, nu = 4, nc = 3, ns = 2, nct = 1;
		DynamicOcpSize const size {g, nx, nu, nc, ns, nct, true};

		for (auto v : vertices(g))
		{
			EXPECT_EQ(size.nx(v), in_degree(v, g) > 0 ? nx : 0) << " at vertex " << v;
			EXPECT_EQ(size.nu(v), out_degree(v, g) > 0 ? nu : 0) << " at vertex " << v;
			EXPECT_EQ(size.nc(v), out_degree(v, g) > 0 ? nc : nct) << " at vertex " << v;
		}
	}


	TEST(DynamicOcpSizeTest, testRobustMpc)
	{
		OcpTree const g(4, 2, 2);

		size_t const nx = 5, nu = 4, nc = 3, ns = 2, nct = 1;
		DynamicOcpSize const size {g, nx, nu, nc, ns, nct, true};

		for (auto v : vertices(g))
		{
			EXPECT_EQ(size.nx(v), in_degree(v, g) > 0 ? nx : 0) << " at vertex " << v;
			EXPECT_EQ(size.nu(v), out_degree(v, g) > 0 ? nu : 0) << " at vertex " << v;
			EXPECT_EQ(size.nc(v), out_degree(v, g) > 0 ? nc : nct) << " at vertex " << v;
		}
	}


	// TEST(DynamicOcpSizeTest, test_ocpSizeNominalMhe)
	// {
	// 	OcpTree const g(4);

	// 	size_t const nx = 5, nw = 4, nc = 3, ns = 2;
	// 	auto const size_map = ocpSizeNominalMhe(num_vertices(g) - 1, nx, nw, nc, ns);

	// 	DynamicOcpSize const size_root(nx, nw, nc, ns);
	// 	DynamicOcpSize const size_leaf(nx, 0, nc, ns);
	// 	DynamicOcpSize const size_other(nx, nw, nc, ns);

	// 	EXPECT_EQ(forcePrint(get(size_map, vertex(0, g))), forcePrint(size_root));

	// 	for (size_t i = 1; i < 3; ++i)
	// 		EXPECT_EQ(forcePrint(get(size_map, vertex(i, g))), forcePrint(size_other)) << i;

	// 	for (size_t i = 3; i < 4; ++i)
	// 		EXPECT_EQ(forcePrint(get(size_map, vertex(i, g))), forcePrint(size_leaf)) << i;
	// }
}