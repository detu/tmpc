#include <tmpc/json/JsonQp.hpp>
#include <tmpc/json/JsonBlaze.hpp>

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <gtest/gtest.h>


namespace tmpc :: testing
{
    TEST(JsonQpTest, test_ocpQpToFromJson)
    {
        using Kernel = BlazeKernel<double>;
        using Stage = OcpQp<Kernel>;

        Stage stage0 {OcpSize {3, 2, 0}, 2};
        Stage stage1 {OcpSize {2, 1, 0}, 0};
        
        randomize(stage0);
        randomize(stage1);

        std::vector<Stage> qp;
        qp.push_back(stage0);
        qp.push_back(stage1);

        json j = ocpQpToJson(qp);

        std::cout << std::setw(4) << j << std::endl;
    }
}
