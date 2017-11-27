#pragma once

#include <tmpc/mpc/MpcOcpSize.hpp>	// for mpcOcpSize()

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include <tmpc/test_tools.hpp>

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>

namespace tmpc :: testing
{
	template <typename WS>
	class SoftConstraintsTest 
	: 	public ::testing::Test
	{
	protected:
		using Workspace = WS;
		using Kernel = typename WS::Kernel;
		using Real = typename Workspace::Real;
		using Vector = DynamicVector<Kernel>;
		using Matrix = DynamicMatrix<Kernel>;
	};


	TYPED_TEST_CASE_P(SoftConstraintsTest);


	TYPED_TEST_P(SoftConstraintsTest, test1dSmallPenalty)
	{
		typename TestFixture::Workspace workspace {OcpSize {1, 0, 0, 1}, OcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.Q(1.);
		stage0.q(-2.);
		stage0.softConstraints({0}, 1e-1, 1e-1);
		stage0.stateBounds(-5., 1.);
		
		workspace.solve();
		auto const solution = workspace.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_PRED2(MatrixApproxEquality(1e-5), solution[0].x(), (Vector {1.81818}));
	}


	TYPED_TEST_P(SoftConstraintsTest, test1dBigPenalty)
	{
		typename TestFixture::Workspace workspace {OcpSize {1, 0, 0, 1}, OcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.Q(1.);
		stage0.q(-2.);
		stage0.softConstraints({0}, 1e+1, 1e+1);
		stage0.stateBounds(-5., 1.);
		
		workspace.solve();
		auto const solution = workspace.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_PRED2(MatrixApproxEquality(1e-5), solution[0].x(), (Vector {1.}));
	}


	TYPED_TEST_P(SoftConstraintsTest, test2stage1dInfeasible)
	{
		typename TestFixture::Workspace workspace {OcpSize {1, 0, 0, 1}, OcpSize {1, 0, 0}};	
		
		auto problem = workspace.problem();
		
		problem[0]
		.Q(1.)
		.q(2.)
		.shootingEquality(1., 0., 0.)
		.softConstraints({0}, 1e+3, 1e+3)
		.stateBounds(-5., -1.);
	
		problem[1]
		.Q(1.)
		.q(-2.)
		.stateBounds(1., 5.);
		
		workspace.solve();
		auto const solution = workspace.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_PRED2(MatrixApproxEquality(1e-5), solution[0].x(), (Vector {1.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-5), solution[1].x(), (Vector {1.}));
	}


	REGISTER_TYPED_TEST_CASE_P(SoftConstraintsTest,
		test1dSmallPenalty,
		test1dBigPenalty,
		test2stage1dInfeasible
	);
}