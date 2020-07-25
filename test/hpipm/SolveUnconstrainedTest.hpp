#pragma once

#include <tmpc/ocp/DynamicOcpKktValue.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/QpSolverTraits.hpp>
#include <tmpc/testing/qp/IsSolution.hpp>

#include <tmpc/Testing.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
    template <typename S>
    class SolveUnconstrainedTest 
    :   public Test
    {
    protected:
        using Solver = S;
		using Real = typename Solver::Real;
		using Qp = DynamicOcpQp<Real>;
		using Solution = DynamicOcpSolution<Real>;
        using Vector = blaze::DynamicVector<Real, blaze::columnVector>;
        using Matrix = blaze::DynamicMatrix<Real>;

		static constexpr bool supportsTreeProblems_ = SupportsTreeProblems<Solver>::value;
    };


    TYPED_TEST_SUITE_P(SolveUnconstrainedTest);

	TYPED_TEST_P(SolveUnconstrainedTest, testNominal)
	{
		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;
		size_t constexpr NT = 2;

		OcpTree const g(NT + 1);
		DynamicOcpSize const size_map {g, NX, NU, NC, 0, NCT, false};
		typename TestFixture::Qp qp {size_map};
		typename TestFixture::Solution sol {size_map};
		typename TestFixture::Solver ws {size_map};
		randomize(qp);
		ws.solveUnconstrained(qp, sol);

		EXPECT_TRUE(isSolution(sol, qp, 1e-9));
	}


	TYPED_TEST_P(SolveUnconstrainedTest, testTree)
	{
		using Solver = typename TestFixture::Solver;

		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;

		OcpTree const g(3, 1, 2);
		DynamicOcpSize const size_map {g, NX, NU, NC, 0, NCT, false};
		typename TestFixture::Qp qp {size_map};
		typename TestFixture::Solution sol {size_map};
		randomize(qp);

		if (this->supportsTreeProblems_)
		{
			Solver solver {size_map};
			solver.solveUnconstrained(qp, sol);
			EXPECT_TRUE(isSolution(sol, qp, 1e-9));
		}
		else
		{
			EXPECT_THROW(Solver {size_map}, std::invalid_argument);
		}
	}


    REGISTER_TYPED_TEST_SUITE_P(SolveUnconstrainedTest
        , testNominal
		, testTree
    );
}