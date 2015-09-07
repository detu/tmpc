#include <MotionPlatformModelPredictiveController.hpp>
#include <MotionPlatformX.hpp>
#include <CyberMotion.hpp>

#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

//int not_main(int argc, char * argv[])
TEST(mpc_test, mpc_test_case)
{
	//const auto platform = std::make_shared<rtmc::MotionPlatformX>();
	const auto platform = std::make_shared<mpmc::CyberMotion>();
	const double Ts = 0.05;
	const unsigned Nt = 10;
	const unsigned simulation_steps = 20;
	const double g = 9.81;
	const double a_max = 1.0;
	const double freq = 1.0;
	std::ofstream out("out.txt");

	mpmc::MotionPlatformModelPredictiveController controller(platform, Ts, Nt);
	controller.setWashoutFactor(0.1);

	Eigen::VectorXd x0(controller.nX());
	x0.fill(0.);
	platform->getDefaultAxesPosition(x0.data());
	controller.InitWorkingPoint(x0);

	Eigen::MatrixXd y_ref(platform->getOutputDim(), Nt);
	y_ref.fill(0.);

	Eigen::VectorXd u(platform->getInputDim());
	u.fill(0.);

	for (unsigned i = 0; i < simulation_steps; ++i)
	{
		std::cout << "step " << i << ", x0=" << x0 << std::endl;

		for (unsigned k = 0; k < y_ref.cols(); ++k)
		{
			y_ref(0, k) = sin(2 * M_PI * freq * i * Ts);
			//y_ref(0, k) = 1;
			y_ref(2, k) = -g;
		}

		controller.setReference(y_ref);
		controller.EmbedInitialValue(x0.data());
		//controller.PrintQP(std::cout);

		try
		{
			controller.Solve();
		}
		catch (const std::runtime_error& e)
		{
			{
				std::ofstream os("failed_qp.m");
				controller.PrintQP_MATLAB(os);
			}

			{
				std::ofstream os("failed_qp.cpp");
				controller.PrintQP_C(os);
			}

			throw e;
		}

		controller.getWorkingU(0, u.data());
		controller.PrepareForNext();

		std::cout << "\tu = " << u << std::endl;
		out << y_ref(0, 0) << "\t" << u.transpose() << std::endl << std::flush;

		auto q = x0.topRows(platform->getNumberOfAxes());
		auto v = x0.bottomRows(platform->getNumberOfAxes());
		x0 << q + v * Ts + u * Ts * Ts / 2., v + u * Ts;
	}
}
