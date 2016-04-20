#include "CyberMotion.hpp"

#include <array>
#include <cassert>
#include <limits>
#include <algorithm>
#include <stdexcept>

#include "cybermotion_generated.h"

namespace mpmc
{
	const unsigned CyberMotion::numberOfAxes;

	CyberMotion::CyberMotion()
	:	MotionPlatform(numberOfAxes)
	,	_ode(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_ode))
	,	_output(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_output))
	{
	}

	void CyberMotion::Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const
	{
		if (x.size() != getStateDim())
			throw std::invalid_argument("CyberMotion::Output(): x has wrong size.");

		if (u.size() != getInputDim())
			throw std::invalid_argument("CyberMotion::Output(): u has wrong size.");

		if (y.size() != getOutputDim())
			throw std::invalid_argument("CyberMotion::Output(): y has wrong size.");

		if (!(C.rows() == getOutputDim() && C.cols() == getStateDim()))
			throw std::invalid_argument("CyberMotion::Output(): C has wrong size.");

		if (!(D.rows() == getOutputDim() && D.cols() == getInputDim()))
			throw std::invalid_argument("CyberMotion::Output(): D has wrong size.");

		if (!(0.5700 <= x(7) && x(7) <= 1.1603))
		{
			std::stringstream msg;
			msg << "CyberMotion::Output(): axis 7 (cabin) position is out of the [0.5700, 1.1603] range. "
					<< "Only this range is currently supported (both rollers on the curve). The value given was " << std::to_string(x(7));
			throw std::out_of_range(msg.str());
		}

		// Call CasADi-generated code.
		_output({x.data(), u.data(), nullptr}, {y.data(), C.data(), D.data()});
	}

	Eigen::VectorXd CyberMotion::getDefaultAxesPosition() const
	{
		// The "agile" position.
		Eigen::VectorXd result(numberOfAxes);
		result << 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991;
		return result;
	}

	void CyberMotion::getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const
	{
		if (q_min.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): q_min has wrong size.");

		if (q_max.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): q_max has wrong size.");

		if (v_min.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): v_min has wrong size.");

		if (v_max.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): v_max has wrong size.");

		if (u_min.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): u_min has wrong size.");

		if (u_max.size() != numberOfAxes)
			throw std::invalid_argument("CyberMotion::getAxesLimits(): u_max has wrong size.");

		static const double data[numberOfAxes][6] = {
			{ 0.2, 9.3, -1.47, 1.47, -1.0780, 1.0780 },
			{ -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), -1.1802, 1.1802, -1.6762, 1.6762 },
			{ -2.1447, -0.8727, -0.9749, 0.9749, -1.1973, 1.1973 },
			{ -0.7540, 1.5415, -1.1802, 1.1802, -2.1893, 2.1893 },
			{ -3.0159, 3.0159, -1.2999, 1.2999, -0.5644, 0.5644 },
			{ -0.8727, 0.8727, -1.2999, 1.2999, -1.6249, 1.6249 },
			{ -3.0159, 3.0159, -2.0525, 2.0525, -1.3170, 1.3170 },
			// Full-range position for cabin axis is { 0.1938, 1.4839, -0.2450, 0.2450, -0.9800, 0.9800 }.
			// Below, the cabin axis position is limited to both rollers on the curve.
			// This is the only case implemented in current CasADi-generated code.
			{ 0.5700, 1.1603, -0.2450, 0.2450, -0.9800, 0.9800 }
		};

		for (int i = 0; i < numberOfAxes; ++i)
		{
			q_min(i) = data[i][0];	q_max(i) = data[i][1];
			v_min(i) = data[i][2];	v_max(i) = data[i][3];
			u_min(i) = data[i][4];	u_max(i) = data[i][5];
		}
	}
}
