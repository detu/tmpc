#include <tmpc/qp/QuadraticProblem.hpp>
#include <tmpc/qp/Printing.hpp>

#include <tmpc/BlazeKernel.hpp>

#include <vector>
#include <iostream>

int main(int, char **)
{
	using namespace tmpc;

	using Kernel = BlazeKernel<double>;
	using Stage = QuadraticProblemStage<Kernel>;

	Stage stage0 {QpSize {3, 2, 0}, 2};
	
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

	Stage stage1 {QpSize {2, 1, 0}, 0};

	stage1.gaussNewtonCostApproximation(
		Kernel::DynamicVector {0.1, 0.2},
		Kernel::DynamicMatrix {
			{0.3, 0.4},
			{0.5, 0.6}
		},
		Kernel::DynamicMatrix {
			{0.7},
			{0.8}
		}
	);

	stage0.linearizedShootingEquality(
		Kernel::DynamicVector {0.9, 1.0},
		Kernel::DynamicMatrix {
			{1.1, 1.2, 1.3},
			{1.4, 1.5, 1.6}
		},
		Kernel::DynamicMatrix {
			{1.7, 1.8},
			{1.9, 2.0}
		},
		Kernel::DynamicVector {2.1, 2.2}
	);

	stage0.relativeBounds(
		Kernel::DynamicVector {0.1, 0.2, 0.3},	// x
		Kernel::DynamicVector {0.4, 0.5},	// u
		Kernel::DynamicVector {-1., -2., -3.},	// lx
		Kernel::DynamicVector {-4., -5.},	// lu
		Kernel::DynamicVector {1., 2., 3.},	// ux
		Kernel::DynamicVector {4., 5.}	// uu
	);

	/*
	blaze::DynamicVector<double> res {0.1, 0.2};
	blaze::DynamicMatrix<double> C {
		{0.3, 0.4},
		{0.5, 0.6}
	};

	blaze::DynamicVector<double> v;
	v = trans(C) * res;
	*/

	std::vector<Stage> qp;
	qp.push_back(stage0);
	qp.push_back(stage1);
	
	std::copy(qp.begin(), qp.end(), std::ostream_iterator<Stage>(std::cout, "\n"));

	return 0;
}