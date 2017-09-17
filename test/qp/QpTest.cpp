#include <tmpc/qp/QuadraticProblem.hpp>

#include <tmpc/BlazeKernel.hpp>

#include "../gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

namespace tmpc :: testing
{
	TEST(QpTest, testStageConstructor)
	{
		using Real = double;
		using Kernel = BlazeKernel<Real>;
		using Stage = QuadraticProblemStage<Kernel>;

		std::vector<Stage> qp;

		Stage s { QpSize {3, 2, 0}, 0 };
		s.Q(Kernel::DynamicMatrix {
			{1., 0.},
			{0., 1.}
		});
	}	
}