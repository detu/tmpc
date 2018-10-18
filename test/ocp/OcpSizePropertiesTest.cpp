#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <tmpc/test_tools.hpp>

#include <vector>


namespace tmpc :: testing
{
	class OcpSizePropertiesTest
	: 	public ::testing::Test
	{
	public:
		OcpSizePropertiesTest()
		{			
		}


	protected:
		size_t const N_ = 5;
		std::vector<size_t> const ns_ = { 2, 1, 1, 0, 0, };
		
		std::vector<OcpSize> const sz_ = 
		{
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {1, 0},
			OcpSize {4, 0},
		};

		OcpGraph const g_ = ocpGraphFromOutDegree(ns_.begin());
		iterator_property_map<std::vector<OcpSize>::const_iterator, decltype(vertexIndex(g_))> const sizeMap_ {sz_.begin(), vertexIndex(g_)};
	};


	TEST_F(OcpSizePropertiesTest, test_size_Q)
	{
		auto const map = size_Q(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nx(), sz_[v].nx()));
	}


	TEST_F(OcpSizePropertiesTest, test_size_R)
	{
		auto const map = size_R(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nu(), sz_[v].nu()));
	}


	TEST_F(OcpSizePropertiesTest, test_size_S)
	{
		auto const map = size_S(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nu(), sz_[v].nx()));
	}


	TEST_F(OcpSizePropertiesTest, test_size_x)
	{
		auto const map = size_x(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), sz_[v].nx());
	}


	TEST_F(OcpSizePropertiesTest, test_size_u)
	{
		auto const map = size_u(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), sz_[v].nu());
	}


	TEST_F(OcpSizePropertiesTest, test_size_C)
	{
		auto const map = size_C(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nc(), sz_[v].nx()));
	}


	TEST_F(OcpSizePropertiesTest, test_size_D)
	{
		auto const map = size_D(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nc(), sz_[v].nu()));
	}


	TEST_F(OcpSizePropertiesTest, test_size_d)
	{
		auto const map = size_d(sizeMap_);
	
		for (auto v : make_iterator_range(vertices(g_)))
			EXPECT_EQ(get(map, v), sz_[v].nc());
	}


	TEST_F(OcpSizePropertiesTest, test_size_A)
	{
		auto const map = size_A(sizeMap_, g_);
	
		for (auto e : make_iterator_range(edges(g_)))
		{
			auto const u = source(e, g_);
			auto const v = target(e, g_);
			EXPECT_EQ(get(map, e), std::pair(sz_[v].nx(), sz_[u].nx()));
		}
	}


	TEST_F(OcpSizePropertiesTest, test_size_B)
	{
		auto const map = size_B(sizeMap_, g_);
	
		for (auto e : make_iterator_range(edges(g_)))
		{
			auto const u = source(e, g_);
			auto const v = target(e, g_);
			EXPECT_EQ(get(map, e), std::pair(sz_[v].nx(), sz_[u].nu()));
		}
	}


	TEST_F(OcpSizePropertiesTest, test_size_b)
	{
		auto const map = size_b(sizeMap_, g_);
	
		for (auto e : make_iterator_range(edges(g_)))
		{
			auto const u = source(e, g_);
			auto const v = target(e, g_);
			EXPECT_EQ(get(map, e), sz_[v].nx());
		}
	}


	TEST(OcpSizeTest, test_ocpSizeNominalMpc)
	{
		OcpGraph const g = ocpGraphLinear(4);

		size_t const nx = 5, nu = 4, nc = 3, ns = 2, nct = 1;
		auto const size_map = ocpSizeNominalMpc(num_vertices(g) - 1, nx, nu, nc, ns, nct, true);

		OcpSize const size_root(0, nu, nc, ns);
		OcpSize const size_leaf(nx, 0, nct, ns);
		OcpSize const size_other(nx, nu, nc, ns);

		EXPECT_EQ(get(size_map, vertex(0, g)), size_root);

		for (size_t i = 1; i < 3; ++i)
			EXPECT_EQ(get(size_map, vertex(i, g)), size_other) << i;

		for (size_t i = 3; i < 4; ++i)
			EXPECT_EQ(get(size_map, vertex(i, g)), size_leaf) << i;
	}


	TEST(OcpSizeTest, test_ocpSizeRobustMpc)
	{
		OcpGraph const g = ocpGraphRobustMpc(4, 2, 2);

		size_t const nx = 5, nu = 4, nc = 3, ns = 2, nct = 1;
		auto const size_map = ocpSizeRobustMpc(g, nx, nu, nc, ns, nct, true);

		OcpSize const size_root(0, nu, nc, ns);
		OcpSize const size_leaf(nx, 0, nct, ns);
		OcpSize const size_other(nx, nu, nc, ns);

		EXPECT_EQ(get(size_map, vertex(0, g)), size_root);

		for (size_t i = 1; i < 7; ++i)
			EXPECT_EQ(get(size_map, vertex(i, g)), size_other) << i;

		for (size_t i = 7; i < 11; ++i)
			EXPECT_EQ(get(size_map, vertex(i, g)), size_leaf) << i;
	}
}