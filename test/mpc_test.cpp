#include <ModelPredictiveController.hpp>
#include <CyberMotionOCP.hpp>
#include <qp/CondensingSolver.hpp>

#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

TEST(mpc_test, mpc_test_case)
{
	typedef mpmc::CyberMotionOCP OCP;
	typedef camels::CondensingSolver<OCP::NX, OCP::NU, OCP::NC, OCP::NCT> QPSolver;
	typedef camels::ModelPredictiveController<OCP, QPSolver> Controller;

	const double Ts = 0.05;
	const unsigned Nt = 1;
	const unsigned simulation_steps = 20;
	const double g = 9.81;
	const double a_max = 1.0;
	const double freq = 1.0;
	std::ofstream out("out.txt");

	OCP ocp(Nt);
	ocp.setWashoutFactor(0.1);

	Controller controller(ocp, Ts);

	auto x0 = ocp.getDefaultState();
	controller.InitWorkingPoint(x0);

	OCP::InputVector u;
	u.fill(0.);

	for (unsigned i = 0; i < simulation_steps; ++i)
	{
		std::cout << "step " << i << ", x0=" << x0 << std::endl;

		for (unsigned k = 0; k < Nt + 1; ++k)
		{
			OCP::OutputVector y_ref = OCP::OutputVector::Zero();
			y_ref[0] = sin(2 * M_PI * freq * i * Ts);
			//y_ref(0, k) = 1;
			y_ref[2] = -g;

			ocp.setReference(k, y_ref);

			if (k == 0)
				out << y_ref[0] << "\t";
		}

		controller.EmbedInitialValue(x0);
		//controller.PrintQP(std::cout);

		try
		{
			controller.Solve();
		}
		catch (const QPSolver::SolveException& e)
		{
			{
				std::ofstream os("failed_qp.m");
				controller.PrintQP_MATLAB(os);
				Print_MATLAB(os, e.getCondensedQP(), "cond_qp");
			}

			{
				std::ofstream os("failed_qp.cpp");
				controller.PrintQP_C(os);
			}

			throw;
		}
		catch (const std::runtime_error&)
		{
			{
				std::ofstream os("failed_qp.m");
				controller.PrintQP_MATLAB(os);
			}

			{
				std::ofstream os("failed_qp.cpp");
				controller.PrintQP_C(os);
			}

			throw;
		}

		u = controller.getWorkingU(0);
		controller.PrepareForNext();

		std::cout << "\tu = " << u << std::endl;
		out << u.transpose() << std::endl << std::flush;

		auto q = x0.topRows<OCP::NU>();
		auto v = x0.bottomRows<OCP::NU>();
		x0 << q + v * Ts + u * Ts * Ts / 2., v + u * Ts;
	}
}
