#include <tmpc/qp/HpmpcWorkspace.hpp>
#include <tmpc/qp/Printing.hpp>

#include <tmpc/BlazeKernel.hpp>
#include <tmpc/Math.hpp>

#include <vector>
#include <iostream>

int main(int, char **)
{
	using namespace tmpc;

	using Kernel = BlazeKernel<double>;
	using Workspace = HpmpcWorkspace<Kernel>;

	Workspace workspace {QpSize {3, 0, 0}, QpSize {0, 0, 0}};
	
	auto& stage0 = workspace.problem()[0];
	stage0.gaussNewtonCostApproximation(
		DynamicVector<double> {1., 2., 42.},
		IdentityMatrix<double> {3u},
		DynamicMatrix<double> {3u, 0u}
	);
	stage0.bounds(-inf<double>(), -inf<double>(), inf<double>(), inf<double>());

	workspace.solve();

	std::cout << workspace.solution()[0].x() << std::endl;

	return 0;
}