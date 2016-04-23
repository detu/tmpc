#include "MotionCueingController.hpp"

#include <Eigen/Dense>

namespace mpmc
{
	MotionCueingController::MotionCueingController(CyberMotionOCP const& ocp, double sample_time, unsigned Nt) :
		camels::ModelPredictiveController<CyberMotionOCP>(CyberMotionOCP::NX, CyberMotionOCP::NU, 2 * CyberMotionOCP::NU, 2 * CyberMotionOCP::NU, sample_time, Nt),
		_ocp(ocp),
		_Nq(CyberMotionOCP::NU),
		_Ny(CyberMotionOCP::NY),
		_yRef(Eigen::Index(CyberMotionOCP::NY), Nt),
		_washoutPosition(Eigen::Index(CyberMotionOCP::NU)),
		_washoutFactor(0.),
		_errorWeight(Eigen::Index(CyberMotionOCP::NY))
	{
		// Initialize limits.
		auto const x_min = _ocp.getStateMin();
		auto const x_max = _ocp.getStateMax();
		auto const u_min = _ocp.getInputMin();
		auto const u_max = _ocp.getInputMax();

		setXMin(x_min);	// TODO: remove boundary constraints for axes position; they must be handled by the non-linear SR-constraints.
		setXMax(x_max);
		setTerminalXMin(x_min);
		setTerminalXMax(x_max);
		setUMin(u_min);
		setUMax(u_max);

		// Initialize weights
		_errorWeight.fill(1.);
	}

	void MotionCueingController::LagrangeTerm(const Eigen::VectorXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) const
	{
		LagrangeTerm(i, z, g, H);
	}
	
	void MotionCueingController::LagrangeTerm(unsigned i, StateInputVector const& z, Eigen::VectorXd& g, Eigen::MatrixXd& H) const
	{
		// Output vector and derivatives.
		OutputVector y;
		OutputJacobianMatrix G;	// G = [C, D]
		_ocp.Output(i, z, y, G);

		const auto W = _errorWeight.cwiseAbs2().asDiagonal();

		// H = G^T W G + \mu I
		H = G.transpose() * W * G;

		/*
		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H.topLeftCorner(nX(), nX()) += _washoutFactor * MatrixXd::Identity(nX(), nX());
		*/

		// g = 2 * (y_bar - y_hat)^T * W * G
		g = (y - _yRef.col(i)).transpose() * W * G;

		/*
		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g.topRows(nX()) += _washoutFactor * (z.topRows(nX()) - getWashoutState());
		*/
	}

	void MotionCueingController::MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) const
	{
		using namespace Eigen;

		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MatrixXd::Identity(nX(), nX()) * getNumberOfIntervals();

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState()) * getNumberOfIntervals();
	}

	void MotionCueingController::setReference(unsigned i, const Eigen::VectorXd& y_ref)
	{
		_yRef.col(i) = y_ref;
	}

	void MotionCueingController::setReference(const Eigen::MatrixXd& y_ref)
	{
		if (!(y_ref.rows() == nY() && y_ref.cols() == getNumberOfIntervals()))
			throw std::invalid_argument("MotionCueingController::setReference(): y_ref has wrong size.");

		_yRef = y_ref;
	}

	Eigen::VectorXd MotionCueingController::getWashoutState() const
	{
		// The "washout" state.
		Eigen::VectorXd x_washout(nX());
		x_washout << _washoutPosition, Eigen::VectorXd::Zero(_Nq);

		return x_washout;
	}

	void MotionCueingController::Integrate(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x_next, Eigen::MatrixXd& A, Eigen::MatrixXd& B) const
	{
		getStateSpaceA(A);
		getStateSpaceB(B);

		x_next = A * x + B * u;		
	}

	void MotionCueingController::getStateSpaceA( Eigen::MatrixXd& A ) const
	{
		using namespace Eigen;

		A << MatrixXd::Identity(_Nq, _Nq), getSampleTime() * MatrixXd::Identity(_Nq, _Nq),
			MatrixXd::Zero(_Nq, _Nq), MatrixXd::Identity(_Nq, _Nq);
	}

	void MotionCueingController::getStateSpaceB( Eigen::MatrixXd& B ) const
	{
		using namespace Eigen;

		B << std::pow(getSampleTime(), 2) / 2. * MatrixXd::Identity(_Nq, _Nq),
			getSampleTime() * MatrixXd::Identity(_Nq, _Nq);
	}

	unsigned MotionCueingController::nY() const
	{
		return _Ny;
	}

	void MotionCueingController::setErrorWeight(const Eigen::VectorXd& val)
	{
		if (val.size() != nY())
			throw std::invalid_argument("MotionCueingController::setErrorWeight(): val has invalid size.");

		_errorWeight = val;
	}

	const Eigen::VectorXd& MotionCueingController::getErrorWeight() const
	{
		return _errorWeight;
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
