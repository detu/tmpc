#include "MotionPlatformModelPredictiveController.hpp"

namespace mpmc
{
	MotionPlatformModelPredictiveController::MotionPlatformModelPredictiveController(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt) :
		camels::MPC_Controller(platform->getStateDim(), platform->getInputDim(), 2 * platform->getNumberOfAxes(), 2 * platform->getNumberOfAxes(), sample_time, Nt),
		_platform(platform),
		_Nq(platform->getNumberOfAxes()),
		_Ny(platform->getOutputDim()),
		_yRef(platform->getOutputDim(), Nt),
		_washoutPosition(platform->getNumberOfAxes()),
		_washoutFactor(0.),
		_errorWeight(platform->getOutputDim())
	{
		// Initialize limits.
		Eigen::VectorXd x_min(nX()), x_max(nX()), u_min(nU()), u_max(nU());
		platform->getAxesLimits(x_min.data(), x_max.data(), x_min.data() + _Nq, x_max.data() + _Nq, u_min.data(), u_max.data());

		setXMin(x_min);	// TODO: remove boundary constraints for axes position; they must be handled by the non-linear SR-constraints.
		setXMax(x_max);
		setTerminalXMin(x_min);
		setTerminalXMax(x_max);
		setUMin(u_min);
		setUMax(u_max);

		// Initialize weights
		_errorWeight.fill(1.);

		// Initialize washout position.
		platform->getDefaultAxesPosition(_washoutPosition.data());
	}

	void MotionPlatformModelPredictiveController::LagrangeTerm(const Eigen::VectorXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) const
	{
		using namespace Eigen;
		
		// Output vector and derivatives.
		// Output() returns derivative matrices in column-major format.
		MatrixXd ssC(_Ny, nX());
		MatrixXd ssD(_Ny, nU());
		VectorXd y(_Ny);
		_platform->Output(z.data(), z.data() + nX(), y.data(), ssC.data(), ssD.data());

		// G = [C, D]
		MatrixXd G(_Ny, nX() + nU());
		G << ssC, ssD;

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

	void MotionPlatformModelPredictiveController::MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) const
	{
		using namespace Eigen;

		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MatrixXd::Identity(nX(), nX()) * getNumberOfIntervals();

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState()) * getNumberOfIntervals();
	}

	void MotionPlatformModelPredictiveController::setReference(unsigned i, const Eigen::VectorXd& y_ref)
	{
		_yRef.col(i) = y_ref;
	}

	void MotionPlatformModelPredictiveController::setReference(const Eigen::MatrixXd& y_ref)
	{
		if (!(y_ref.rows() == nY() && y_ref.cols() == getNumberOfIntervals()))
			throw std::invalid_argument("MotionPlatformModelPredictiveController::setReference(): y_ref has wrong size.");

		_yRef = y_ref;
	}

	Eigen::VectorXd MotionPlatformModelPredictiveController::getWashoutState() const
	{
		// The "washout" state.
		Eigen::VectorXd x_washout(nX());
		x_washout << _washoutPosition, Eigen::VectorXd::Zero(_Nq);

		return x_washout;
	}

	void MotionPlatformModelPredictiveController::Integrate(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x_next, Eigen::MatrixXd& A, Eigen::MatrixXd& B) const
	{
		getStateSpaceA(A);
		getStateSpaceB(B);

		x_next = A * x + B * u;		
	}

	void MotionPlatformModelPredictiveController::getStateSpaceA( Eigen::MatrixXd& A ) const
	{
		using namespace Eigen;

		A << MatrixXd::Identity(_Nq, _Nq), getSampleTime() * MatrixXd::Identity(_Nq, _Nq),
			MatrixXd::Zero(_Nq, _Nq), MatrixXd::Identity(_Nq, _Nq);
	}

	void MotionPlatformModelPredictiveController::getStateSpaceB( Eigen::MatrixXd& B ) const
	{
		using namespace Eigen;

		B << std::pow(getSampleTime(), 2) / 2. * MatrixXd::Identity(_Nq, _Nq),
			getSampleTime() * MatrixXd::Identity(_Nq, _Nq);
	}

	unsigned MotionPlatformModelPredictiveController::nY() const
	{
		return _Ny;
	}

	void MotionPlatformModelPredictiveController::setErrorWeight(const Eigen::VectorXd& val)
	{
		if (val.size() != nY())
			throw std::invalid_argument("MotionPlatformModelPredictiveController::setErrorWeight(): val has invalid size.");

		_errorWeight = val;
	}

	const Eigen::VectorXd& MotionPlatformModelPredictiveController::getErrorWeight() const
	{
		return _errorWeight;
	}

	void MotionPlatformModelPredictiveController::setWashoutPosition(const Eigen::VectorXd& val)
	{
		if (val.size() != _Nq)
			throw std::invalid_argument("MotionPlatformModelPredictiveController::setWashoutPosition(): val has invalid size.");

		_washoutPosition = val;
	}

	template<class Vector1, class Matrix, class Vector2>
	void MotionPlatformModelPredictiveController::SRConstraints(const Eigen::MatrixBase<Vector1>& x, Eigen::MatrixBase<Matrix>& D, Eigen::MatrixBase<Vector2>& d_min, Eigen::MatrixBase<Vector2>& d_max) const
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

	void MotionPlatformModelPredictiveController::PathConstraints(unsigned i, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const
	{
		auto Dx = D.leftCols(nX());
		SRConstraints(x, Dx, d_min, d_max);
		D.rightCols(nU()).setConstant(0.);
	}

	void MotionPlatformModelPredictiveController::TerminalConstraints(const Eigen::VectorXd& x, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const
	{
		SRConstraints(x, D, d_min, d_max);
	}
}
