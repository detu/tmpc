#include <tmpc/qp/QpOasesWorkspace.hpp>
#include <tmpc/qp/Printing.hpp>

#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>
#include <tmpc/Math.hpp>

#include <vector>
#include <iostream>

int main(int, char **)
{
	using namespace tmpc;

	using Kernel = BlazeKernel<double>;
	using Workspace = QpOasesWorkspace<Kernel>;

	Workspace workspace {QpSize {3, 0, 0}, QpSize {0, 0, 0}};
	
	auto& stage0 = workspace.problem()[0];
	stage0.gaussNewtonCostApproximation(
		DynamicVector<Kernel> {1., 2., 42.},
		IdentityMatrix<Kernel> {3u},
		DynamicMatrix<Kernel> {3u, 0u}
	);
	stage0.bounds(-inf<double>(), -inf<double>(), inf<double>(), inf<double>());

	for (auto const& s : workspace.problem())
		std::cout << s << std::endl;

	workspace.solve();

	std::cout << workspace.solution()[0].x() << std::endl;

	return 0;
}