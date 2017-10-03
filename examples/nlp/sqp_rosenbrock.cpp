/**
 * \brief Shows how to use tmpc to solve an NLP using SQP method. 
 * 
 * The sample NLP is taken from here: 
 * https://www.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html
 */


#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpPoint.hpp>
#include <tmpc/qp/QpOasesWorkspace.hpp>
#include <tmpc/print/qp/OcpQp.hpp>
#include <tmpc/print/ocp/OcpSolution.hpp>

#include <tmpc/BlazeKernel.hpp>

#include <iostream>
#include <cmath>

using namespace tmpc;

using Real = double;
using Kernel = BlazeKernel<Real>;


/// \brief Defines the objective function, its gradient and Hessian
auto rosenbrock(StaticVector<Kernel, 2> const& x)
{
	return std::make_tuple(
		100. * pow(x[1] - pow(x[0], 2), 2) + pow(1. - x[0], 2),	// objective
		StaticVector<Kernel, 2>	// gradient
		{
			2. * (-1. + x[0] + 200 * pow(x[0], 3) - 200. * x[0] * x[1]), 
			200. * (-pow(x[0], 2) + x[1])
		},
		StaticMatrix<Kernel, 2, 2>	// Hessian
		{
			{2. + 1200. * pow(x[0], 2) - 400. * x[1], -400. * x[0]}, 
			{-400. * x[0], 200.}
		}
	);
}


/// \brief Defines the constraint function, its gradient and Hessian
auto unitdisk(StaticVector<Kernel, 2> const& x)
{
	return std::make_tuple(
		pow(x[0], 2) + pow(x[1], 2),	// the constraint is x[0]^2 + x[1]^2 <= 1
		StaticMatrix<Kernel, 1, 2> {{2. * x[0], 2. * x[1]}},	// gradient
		evaluate(2. * IdentityMatrix<Kernel>(2))	// Hessian
	);
}


int main(int, char **)
{	
	using QpWorkspace = QpOasesWorkspace<Kernel>;
	using Point = OcpPoint<Kernel>;
	
	// Create QP workspace, specify problem size.
	QpWorkspace ws {OcpSize {2, 0, 1}};

	// Set initial point
	std::array<Point, 1> x {Point {{0., 0.}, {}}};

	// Get QP iterator range
	auto qp = ws.problem();

	// Set bounds
	qp[0].lbx(-inf<Real>());
	qp[0].ubx(inf<Real>());
	qp[0].lbd(-inf<Real>());
	qp[0].ubd(1.);

	// Lagrange multiplier
	StaticVector<Kernel, 1> lambda {0.};

	for (int count = 0; count < 10; ++count)
	{
		// Compute quadratic approximation
		auto const f = rosenbrock(x[0].x());
		auto const g = unitdisk(x[0].x());

		// Should it be + or - lambda?
		qp[0].Q(std::get<2>(f) + lambda[0] * std::get<2>(g));
		qp[0].q(std::get<1>(f));
		qp[0].linearizedGeneralConstraints(std::get<0>(g), std::get<1>(g), StaticMatrix<Kernel, 1, 0> {}, -inf<Real>(), 1.);

		std::cout << "x = " << trans(x[0].x()) << ", f = " << std::get<0>(f)
			<< ", g = " << std::get<0>(g) << ", lambda = " << lambda << std::endl;
		//std::cout << "qp[0] = " << std::endl << qp[0] << std::endl;

		// Calculate Newton step
		ws.solve();
		auto const qp_solution = ws.solution();

		//std::cout << "Newton step: " << std::endl << qp_solution[0] << std::endl;

		x[0].x(x[0].x() + qp_solution[0].x());
		lambda = qp_solution[0].lam_ubd();
	}
	

	return 0;
}