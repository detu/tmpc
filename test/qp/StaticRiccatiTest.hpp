/// This test suite is not finished and is excluded from build.

#pragma once

#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>
#include <tmpc/Testing.hpp>
#include <tmpc/testing/qp/IsSolution.hpp>

#include <blaze/Math.h>

namespace tmpc :: testing
{
    template <typename Riccati_>
    class StaticRiccatiTest
    :   public Test
    {
    protected:
		using Riccati = Riccati_;
        using Real = typename Riccati::Real;
    };


    TYPED_TEST_SUITE_P(StaticRiccatiTest);


    TYPED_TEST_P(StaticRiccatiTest, testKktValue)
	{
		using Real = typename TestFixture::Real;

		OcpTree const g {
			3,
			2, 2, 2, 
			1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0
		};
		typename TestFixture::Riccati riccati {g};
		auto const& sz = riccati.size();

		DynamicOcpQp<Real> qp {sz};
		DynamicOcpSolution<Real> sol {sz};
		
		randomize(qp);
		riccati(qp, sol);

		// Set Lagrange multipliers for inequalities to 0
		for (auto u : vertices(g))
		{
			reset(sol.lam_lx(u));
			reset(sol.lam_ux(u));
			reset(sol.lam_lu(u));
			reset(sol.lam_uu(u));
			reset(sol.lam_ld(u));
			reset(sol.lam_ud(u));
		}

		EXPECT_TRUE(isSolution(sol, qp, 1e-10));
	}


    REGISTER_TYPED_TEST_SUITE_P(StaticRiccatiTest
        , testKktValue
    );
}