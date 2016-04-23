#include "MotionCueingController.hpp"

#include <Eigen/Dense>

namespace mpmc
{
	MotionCueingController::MotionCueingController(CyberMotionOCP const& ocp, double sample_time, unsigned Nt) :
		camels::ModelPredictiveController<CyberMotionOCP>(ocp, CyberMotionOCP::NX, CyberMotionOCP::NU, 2 * CyberMotionOCP::NU, 2 * CyberMotionOCP::NU, sample_time, Nt),
		_Nq(CyberMotionOCP::NU),
		_Ny(CyberMotionOCP::NY),
		_washoutPosition(Eigen::Index(CyberMotionOCP::NU)),
		_washoutFactor(0.)
	{
	}
	
	void MotionCueingController::MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) const
	{
		using namespace Eigen;

		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MatrixXd::Identity(nX(), nX()) * getNumberOfIntervals();

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState()) * getNumberOfIntervals();
	}

	Eigen::VectorXd MotionCueingController::getWashoutState() const
	{
		// The "washout" state.
		Eigen::VectorXd x_washout(nX());
		x_washout << _washoutPosition, Eigen::VectorXd::Zero(_Nq);

		return x_washout;
	}

	unsigned MotionCueingController::nY() const
	{
		return _Ny;
	}

	void MotionCueingController::setWashoutPosition(const Eigen::VectorXd& val)
	{
		if (val.size() != _Nq)
			throw std::invalid_argument("MotionCueingController::setWashoutPosition(): val has invalid size.");

		_washoutPosition = val;
	}

	template<class Vector1, class Matrix, class Vector2>
	void MotionCueingController::SRConstraints(const Eigen::MatrixBase<Vector1>& x, Eigen::MatrixBase<Matrix>& D, Eigen::MatrixBase<Vector2>& d_min, Eigen::MatrixBase<Vector2>& d_max) const
	{
		using Eigen::VectorXd;
		using Eigen::MatrixXd;

		const auto q = x.topRows(_Nq);
		const auto v = x.bottomRows(_Nq);
		const VectorXd q_min = getXMin().topRows(_Nq);
		const VectorXd q_max = getXMax().topRows(_Nq);

		/*
		Stoppability-reachability equations:
		v^2 + 2 * (q_max - q) * a_min <= 0
		v^2 + 2 * (q_min - q) * a_max <= 0

		Linearization:
		-inf <= 2 * (-a_min * dq + v * dv) <= -v^2 - 2 * (q_max - q) * a_min
		-inf <= 2 * (-a_min * dq + v * dv) <= -v^2 - 2 * (q_min - q) * a_max
		*/
		D << -2. * MatrixXd(getUMin().asDiagonal()), 2. * MatrixXd(q.asDiagonal()),
			 -2. * MatrixXd(getUMax().asDiagonal()), 2. * MatrixXd(q.asDiagonal());

		const auto inf = VectorXd::Constant(_Nq, std::numeric_limits<double>::infinity());
		d_min << -inf, -inf;
		d_max << -v.cwiseAbs2() - 2. * (q_max - q).cwiseProduct(getUMin()), -v.cwiseAbs2() - 2. * (q_min - q).cwiseProduct(getUMax());
	}

	void MotionCueingController::PathConstraints(unsigned i, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const
	{
		auto Dx = D.leftCols(nX());
		SRConstraints(x, Dx, d_min, d_max);
		D.rightCols(nU()).setConstant(0.);
	}

	void MotionCueingController::TerminalConstraints(const Eigen::VectorXd& x, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const
	{
		SRConstraints(x, D, d_min, d_max);
	}
}
