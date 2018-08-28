#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/core/IteratorRange.hpp>
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
			EXPECT_EQ(get(map, v), std::pair(sz_[v].nx(), sz_[v].nu()));
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
}