#include <MPC_Controller.hpp>
#include <MotionPlatformX.hpp>
#include <CyberMotion.hpp>

#include <iostream>
#include <fstream>

int main(int argc, char * argv[])
{
	//const auto platform = std::make_shared<rtmc::MotionPlatformX>();
	const auto platform = std::make_shared<rtmc::CyberMotion>();
	const double Ts = 0.05;
	const unsigned Nt = 1;
	const unsigned simulation_steps = 200;
	const double g = 9.81;
	const double a_max = 1.0;
	const double freq = 1.0;
	std::ofstream out("out.txt");

	/** Set qpDUNES options */
	qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();
	//qpOptions.maxIter = 100;
	//qpOptions.printLevel = 2;
	//qpOptions.stationarityTolerance = 1.e-6;
	qpOptions.lsType = QPDUNES_LS_BACKTRACKING_LS;

	try
	{
		std::ofstream call_log("qpDUNES_call.log");
		rtmc::MPC_Controller controller(platform, Ts, Nt, &qpOptions, &call_log);
		//controller.qpDUNES_call_log(&call_log);
		controller.InitWorkingPoint();

		Eigen::VectorXd x0(platform->getStateDim());
		x0.fill(0.);
		platform->getDefaultAxesPosition(x0.data());

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

			controller.SetReference(y_ref.data());
			controller.EmbedInitialValue(x0.data());
			//controller.PrintQP(std::cout);

			try
			{
				controller.Solve();
			}
			catch (const std::runtime_error& e)
			{
				controller.PrintQP_MATLAB(std::ofstream("failed_qp.m"));
				controller.PrintQP_C(std::ofstream("failed_qp.cpp"));
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
	catch (const std::runtime_error& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}