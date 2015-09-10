#include <CyberMotion1D.hpp>

#include <array>
#include <cassert>
#include <limits>
#include <algorithm>

namespace mpmc
{
	CyberMotion1D::CyberMotion1D() : 
		MotionPlatform(1)
	{
		_fullState.setConstant(std::numeric_limits<double>::quiet_NaN());
	}

	void CyberMotion1D::Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const
	{
		FullStateVector X = _fullState;
		X(0) = x(0);
		X(CyberMotion::numberOfAxes) = x(1);

		FullInputVector U = FullInputVector::Zero();
		U(0) = u[0];

		Eigen::MatrixXd tmpC(getOutputDim(), _fullCMS.getStateDim());
		Eigen::MatrixXd tmpD(getOutputDim(), _fullCMS.getInputDim());
		_fullCMS.Output(X, U, y, tmpC, tmpD);

		C << tmpC.col(0), tmpC.col(CyberMotion::numberOfAxes);
		D << tmpD.col(0);
	}

	Eigen::VectorXd CyberMotion1D::getDefaultAxesPosition() const
	{
		return _fullCMS.getDefaultAxesPosition().row(0);
	}

	void CyberMotion1D::getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const
	{
		Eigen::VectorXd full_q_min(CyberMotion::numberOfAxes), full_q_max(CyberMotion::numberOfAxes), 
			full_v_min(CyberMotion::numberOfAxes), full_v_max(CyberMotion::numberOfAxes), 
			full_u_min(CyberMotion::numberOfAxes), full_u_max(CyberMotion::numberOfAxes);
		_fullCMS.getAxesLimits(full_q_min, full_q_max, full_v_min, full_v_max, full_u_min, full_u_max);
		
		q_min << full_q_min(0);
		q_max << full_q_max(0);
		v_min << full_v_min(0);
		v_max << full_v_max(0);
		u_min << full_u_min(0);
		u_max << full_u_max(0);
	}

	const CyberMotion1D::FullStateVector& CyberMotion1D::getFullState() const
	{
		return _fullState;
	}

	void CyberMotion1D::setFullState(const FullStateVector& val)
	{
		_fullState = val;
	}

}
