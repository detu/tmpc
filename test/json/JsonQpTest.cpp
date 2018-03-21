#include <tmpc/json/JsonQp.hpp>
#include <tmpc/json/JsonBlaze.hpp>

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <gtest/gtest.h>

#include <tmpc/BlazeKernel.hpp>


namespace tmpc :: testing
{
    TEST(JsonQpTest, test_OcpQpToFromJson)
    {
        using Kernel = BlazeKernel<double>;
        using Stage = OcpQp<Kernel>;

        Stage stage0 {OcpSize {3, 2, 0}, 2};
        
        stage0
        .Q({{1., 0., 0.},
            {0., 2., 0.},
            {0., 0., 3.}})
        .R({{5., 0.},
            {0., 6.}})
        .S({{7., 8.},
            {9., 10.},
            {11., 12.}})
        .q({13., 14., 15.})
        .r({16., 17.});

        json j = stage0;
        std::cout << std::setw(4) << j << std::endl;

        auto stage1 = j.get<Stage>();
    }


    TEST(JsonQpTest, test_toJson)
    {
        using Kernel = BlazeKernel<double>;
        using Stage = OcpQp<Kernel>;

        Stage stage0 {OcpSize {3, 2, 0}, 2};
        
        stage0
        .Q({{1., 0., 0.},
            {0., 2., 0.},
            {0., 0., 3.}})
        .R({{5., 0.},
            {0., 6.}})
        .S({{7., 8.},
            {9., 10.},
            {11., 12.}})
        .q({13., 14., 15.})
        .r({16., 17.});

        Stage stage1 {OcpSize {2, 1, 0}, 0};

        stage1.gaussNewtonCostApproximation(
            DynamicVector<Kernel> {0.1, 0.2},
            DynamicMatrix<Kernel> {
                {0.3, 0.4},
                {0.5, 0.6}
            },
            DynamicMatrix<Kernel> {
                {0.7},
                {0.8}
            }
        );

        stage0.linearizedShootingEquality(
            DynamicVector<Kernel> {0.9, 1.0},
            DynamicMatrix<Kernel> {
                {1.1, 1.2, 1.3},
                {1.4, 1.5, 1.6}
            },
            DynamicMatrix<Kernel> {
                {1.7, 1.8},
                {1.9, 2.0}
            },
            DynamicVector<Kernel> {2.1, 2.2}
        );

        stage0.relativeStateBounds(
            DynamicVector<Kernel> {0.1, 0.2, 0.3},	// x
            DynamicVector<Kernel> {-1., -2., -3.},	// lx
            DynamicVector<Kernel> {1., 2., 3.}		// ux
        );

        stage0.relativeInputBounds(
            DynamicVector<Kernel> {0.4, 0.5},	// u
            DynamicVector<Kernel> {-4., -5.},	// lu
            DynamicVector<Kernel> {4., 5.}	// uu
        );

        std::vector<Stage> qp;
        qp.push_back(stage0);
        qp.push_back(stage1);

        json j = stage0;
    }
}
