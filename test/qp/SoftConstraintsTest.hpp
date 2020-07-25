#pragma once

#include <tmpc/mpc/MpcOcpSize.hpp>	// for mpcOcpSize()
#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Testing.hpp>

#include <iostream>
#include <fstream>


namespace tmpc :: testing
{
	template <typename WS>
	class SoftConstraintsTest 
	: 	public Test
	{
	protected:
		using Workspace = WS;
		using Kernel = typename WS::Kernel;
		using Real = typename Workspace::Real;
		using Vector = DynamicVector<Kernel>;
		using Matrix = DynamicMatrix<Kernel>;
	};


	TYPED_TEST_SUITE_P(SoftConstraintsTest);


	TYPED_TEST_P(SoftConstraintsTest, test1dSmallPenalty)
	{
		typename TestFixture::Workspace workspace {DynamicOcpSize {1, 0, 0, 1}, DynamicOcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.Q(1.);
		stage0.q(-2.);
		stage0.softConstraints({0}, 1e-1, 1e-1);
		stage0.stateBounds(-5., 1.);
		
		workspace.solve();
		auto const solution = workspace.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_TRUE(approxEqual(solution[0].x(), (Vector {1.81818}), 1e-5));
	}


	TYPED_TEST_P(SoftConstraintsTest, test1dBigPenalty)
	{
		typename TestFixture::Workspace workspace {DynamicOcpSize {1, 0, 0, 1}, DynamicOcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.Q(1.);
		stage0.q(-2.);
		stage0.softConstraints({0}, 1e+1, 1e+1);
		stage0.stateBounds(-5., 1.);
		
		workspace.solve();
		auto const solution = workspace.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_TRUE(approxEqual(solution[0].x(), (Vector {1.}), 1e-5));
	}


	TYPED_TEST_P(SoftConstraintsTest, test2stage1dInfeasible)
	{
		typename TestFixture::Workspace workspace {DynamicOcpSize {1, 0, 0, 1}, DynamicOcpSize {1, 0, 0}};	
		
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

		EXPECT_TRUE(approxEqual(solution[0].x(), (Vector {1.}), 1e-5));
		EXPECT_TRUE(approxEqual(solution[1].x(), (Vector {1.}), 1e-5));
	}


	REGISTER_TYPED_TEST_SUITE_P(SoftConstraintsTest,
		test1dSmallPenalty,
		test1dBigPenalty,
		test2stage1dInfeasible
	);
}