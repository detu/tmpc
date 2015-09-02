#include "MotionPlatformModelPredictiveController.hpp"

namespace mpmc
{
	MotionPlatformModelPredictiveController::MotionPlatformModelPredictiveController(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt) :
		MPC_Controller(platform->getStateDim(), platform->getInputDim(), sample_time, Nt),
		_platform(platform),
		_Nq(platform->getNumberOfAxes()),
		_Ny(platform->getOutputDim()),
		_yRef(platform->getOutputDim(), Nt),
		_washoutPosition(platform->getNumberOfAxes()),
		_washoutFactor(0.)
	{
		// Allocate arrays.
		_W.resize(_Ny * _Ny * Nt);
		_G.resize(_Ny * nZ() * Nt);

		// Initialize limits.
		Eigen::VectorXd x_min(nX()), x_max(nX()), u_min(nU()), u_max(nU());
		platform->getAxesLimits(x_min.data(), x_max.data(), x_min.data() + _Nq, x_max.data() + _Nq, u_min.data(), u_max.data());

		setXMin(x_min);
		setXMax(x_max);
		setUMin(u_min);
		setUMax(u_max);

		// Initialize weight matrices
		for (unsigned i = 0; i < Nt; ++i)
			W(i).setIdentity();

		// Initialize washout position.
		platform->getDefaultAxesPosition(_washoutPosition.data());
	}

	void MotionPlatformModelPredictiveController::LagrangeTerm(const Eigen::MatrixXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g)
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

		// H = G^T W G + \mu I
		H = G.transpose() * W(i) * G;

		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H.topLeftCorner(nX(), nX()) += _washoutFactor * MatrixXd::Identity(nX(), nX());

		// g = 2 * (y_bar - y_hat)^T * W * G
		g = (y - _yRef.col(i)).transpose() * W(i) * G;

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g.topRows(nX()) += _washoutFactor * (z.topRows(nX()) - getWashoutState());
	}

	void MotionPlatformModelPredictiveController::MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g)
	{
		using namespace Eigen;

		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MatrixXd::Identity(nX(), nX());

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState());
	}

	MPC_Controller::RowMajorMatrixMap MotionPlatformModelPredictiveController::G(unsigned i)
	{
		assert(i < getNumberOfIntervals());
		return RowMajorMatrixMap(_G.data() + i * _Ny * nZ(), _Ny, nZ());
	}

	void MotionPlatformModelPredictiveController::setReference(unsigned i, const Eigen::VectorXd& y_ref)
	{
		_yRef.col(i) = y_ref;
	}

	void MotionPlatformModelPredictiveController::setReference(const Eigen::MatrixXd& y_ref)
	{
		_yRef = y_ref;
	}

	Eigen::VectorXd MotionPlatformModelPredictiveController::getWashoutState() const
	{
		// The "washout" state.
		Eigen::VectorXd x_washout(nX());
		x_washout << _washoutPosition, Eigen::VectorXd::Zero(_Nq);

		return x_washout;
	}

	MPC_Controller::RowMajorMatrixMap MotionPlatformModelPredictiveController::W(unsigned i)
	{
		assert(i < getNumberOfIntervals());
		return RowMajorMatrixMap(_W.data() + i * _Ny * _Ny, _Ny, _Ny);
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
}