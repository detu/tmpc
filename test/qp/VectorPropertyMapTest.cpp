#include <tmpc/qp/detail/VectorPropertyMap.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <blaze/Blaze.h>

#include <tmpc/test_tools.hpp>
#include <gtest/gtest.h>

#include <vector>


namespace tmpc :: testing
{
    class VectorPropertyMapTest 
    :   public ::testing::Test
    {
    protected:
        VectorPropertyMapTest()
        {
        }


        OcpGraph const graph_ = ocpGraphLinear(3);

        std::vector<OcpSize> const size_ = {
            OcpSize {2, 1, 1},
            OcpSize {5, 4, 2},
            OcpSize {8, 7, 3}
        };
    };


    TEST_F(VectorPropertyMapTest, test_get)
    {
        using Vector = blaze::DynamicVector<double>;
        Vector Q;

        auto const index_map = vertexIndex(graph_);
        auto const size_map = iterator_property_map(size_.begin(), index_map);

        auto const setter = [] (double const * src, int index)
        {

        };

        int expected_index = -1;

        auto const getter = [this, &expected_index] (double * dst, int index)
        {
            EXPECT_EQ(index, expected_index);

            for (size_t i = 0; i < size_[index].nx(); ++i)
                dst[i] = index + i * 0.1;
        };

        auto map = tmpc::detail::makeVectorPropertyMap<OcpVertexDescriptor, Vector>(index_map, size_x(size_map), setter, getter);
        
        expected_index = 0;
        EXPECT_EQ(forcePrint(get(map, 0)), (Vector {0.0, 0.1}));

        expected_index = 1;
        EXPECT_EQ(forcePrint(get(map, 1)), (Vector {1.0, 1.1, 1.2, 1.3, 1.4}));
    }
}
