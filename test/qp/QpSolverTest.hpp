#pragma once

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/QpSolverTraits.hpp>

#include <tmpc/Testing.hpp>
#include <tmpc/testing/qp/IsSolution.hpp>

// #include <tmpc/hpipm/NominalSolver.hpp>
// #include <tmpc/print/ocp/OcpSolution.hpp>
// #include <iostream>


namespace tmpc :: testing
{
	template <typename Solver_>
	class QpSolverTest 
	: 	public Test
	{
	protected:
		using Solver = Solver_;
		using Real = typename Solver::Real;
		using Qp = DynamicOcpQp<Real>;
		using Solution = DynamicOcpSolution<Real>;
		using Vector = blaze::DynamicVector<Real, blaze::columnVector>;
		using Matrix = blaze::DynamicMatrix<Real>;

		static constexpr bool supportsTreeProblems_ = SupportsTreeProblems<Solver>::value;
	};


	TYPED_TEST_SUITE_P(QpSolverTest);

	
	TYPED_TEST_P(QpSolverTest, testSingleVertex)
	{
		size_t constexpr NX = 2;

		DynamicOcpSize const size {OcpTree {0}, {{NX, 0, 0}}};
		typename TestFixture::Qp qp {size};
		typename TestFixture::Solution sol {size};
		typename TestFixture::Solver solver {size};
		
		randomize(qp);
		solver(qp, sol);

		EXPECT_TRUE(isSolution(sol, qp, 1e-9));
	}


	TYPED_TEST_P(QpSolverTest, testSingleVertexGeneralConstraint)
	{
		size_t constexpr NX = 2;
		size_t constexpr NC = 1;

		DynamicOcpSize const size {OcpTree {0}, {{NX, 0, NC}}};
		typename TestFixture::Qp qp {size};
		typename TestFixture::Solution sol {size};
		typename TestFixture::Solver solver {size};
		
		randomize(qp);
		solver(qp, sol);

		EXPECT_TRUE(isSolution(sol, qp, 1e-8));
	}


	TYPED_TEST_P(QpSolverTest, testTwoVertices)
	{
		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;

		DynamicOcpSize const size {OcpTree {1, 0}, 
			{{NX, NU, NC}, {NX, 0, NC}}};
		typename TestFixture::Qp qp {size};
		typename TestFixture::Solution sol {size};
		typename TestFixture::Solver solver {size};
		
		randomize(qp);
		solver(qp, sol);

		EXPECT_TRUE(isSolution(sol, qp, 1e-9));

		// hpipm::NominalSolver<typename TestFixture::Real> nominal_solver {size};
		// typename TestFixture::Solution ref_sol {size};
		// nominal_solver(qp, ref_sol);
		
		// // std::cout << "********      Solution      ******** \n" << sol;
		// // std::cout << "******** Reference solution ******** \n" << ref_sol;
	}

	
	TYPED_TEST_P(QpSolverTest, testSolveNominal)
	{
		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;
		size_t constexpr NT = 2;

		DynamicOcpSize const size (NT, NX, NU, NC, 0, NCT, false);
		typename TestFixture::Qp qp {size};
		typename TestFixture::Solution sol {size};
		typename TestFixture::Solver solver {size};
		
		randomize(qp);
		solver(qp, sol);

		EXPECT_TRUE(isSolution(sol, qp, 1e-9));

		// hpipm::NominalSolver<typename TestFixture::Real> nominal_solver {size};
		// typename TestFixture::Solution ref_sol {size};
		// nominal_solver(qp, ref_sol);

		// std::cout << "********      Solution      ******** \n" << sol;
		// std::cout << "******** Reference solution ******** \n" << ref_sol;
	}


	TYPED_TEST_P(QpSolverTest, testSolveTree)
	{
		using Solver = typename TestFixture::Solver;

		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;
		size_t constexpr NT = 2;

		OcpTree const g {
			3,
			2, 2, 2,
			1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0
		};
		DynamicOcpSize const size {g, NX, NU, NC, 0, NCT, false};
		typename TestFixture::Qp qp {size};
		typename TestFixture::Solution sol {size};
		
		randomize(qp);

		if (this->supportsTreeProblems_)
		{
			Solver solver {size};
			solver(qp, sol);
			EXPECT_TRUE(isSolution(sol, qp, 1e-11));
		}
		else
		{
			EXPECT_THROW(Solver {size}, std::invalid_argument);
		}
	}


	REGISTER_TYPED_TEST_SUITE_P(QpSolverTest
		, testSingleVertex
		, testSingleVertexGeneralConstraint
		, testTwoVertices
		, testSolveNominal
		, testSolveTree
	);
}