#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/print/qp/OcpQp.hpp>
#include <tmpc/print/qp/OcpSolution.hpp>

#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>
#include <tmpc/Math.hpp>

#include <vector>
#include <iostream>

int main(int, char **)
{
	using namespace tmpc;

	using Kernel = BlazeKernel<double>;
	using Workspace = HpipmWorkspace<Kernel>;

	Workspace workspace {OcpSize {1, 0, 0, 1}, OcpSize {1, 0, 0}};	
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
    
    std::cout << "==== Problem ====" << std::endl;
	for (auto const& s : workspace.problem())
		std::cout << s << std::endl;

    workspace.solve();
    auto const solution = workspace.solution();

    std::cout << "==== Solution ====" << std::endl;
    for (int i = 0; i < problem.size(); ++i)
        std::cout << "STAGE " << i << ":" << std::endl << solution[i] << std::endl;

	return 0;
}