#include <tmpc/property_map/VectorPropertyMap.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Blaze.h>

#include <vector>


namespace tmpc :: testing
{
    class VectorPropertyMapTest 
    :   public Test
    {
    protected:
        VectorPropertyMapTest()
        {
        }


        OcpTree const graph_(3);

        std::vector<DynamicOcpSize> const size_ = {
            DynamicOcpSize {2, 1, 1},
            DynamicOcpSize {5, 4, 2},
            DynamicOcpSize {8, 7, 3}
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

        auto map = makeVectorPropertyMap<OcpVertex, Vector>(index_map, size_x(size_map), setter, getter);
        
        expected_index = 0;
        EXPECT_EQ(forcePrint(get(map, 0)), (Vector {0.0, 0.1}));

        expected_index = 1;
        EXPECT_EQ(forcePrint(get(map, 1)), (Vector {1.0, 1.1, 1.2, 1.3, 1.4}));
    }
}
