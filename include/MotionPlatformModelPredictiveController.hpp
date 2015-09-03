#pragma once

#include "MPC_Controller.hpp"
#include "MotionPlatform.hpp"

namespace mpmc
{
	class MotionPlatformModelPredictiveController : public camels::MPC_Controller
	{
	public:
		MotionPlatformModelPredictiveController(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt);

		unsigned nY() const;

		const Eigen::VectorXd& getErrorWeight() const;
		void setErrorWeight(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getWashoutPosition() const { return _washoutPosition; }
		void setWashoutPosition(const Eigen::VectorXd& val);

		double getWashoutFactor() const { return _washoutFactor; }
		void setWashoutFactor(double val) { _washoutFactor = val; }

		void setReference(const Eigen::MatrixXd& py_ref);
		void setReference(unsigned i, const Eigen::VectorXd& py_ref);

	protected:
		void LagrangeTerm(const Eigen::MatrixXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) override;
		void MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) override;
		void Integrate(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x_next, Eigen::MatrixXd& A, Eigen::MatrixXd& B) const override;

	private:
		const std::shared_ptr<MotionPlatform> _platform;
		const unsigned _Nq;
		const unsigned _Ny;

		// Reference trajectory.
		// _yRef stores _Nt vectors of size _Ny.
		// Important: _yRef is column-major.
		Eigen::MatrixXd _yRef;

		// Tracking error weights for all components of output.
		Eigen::VectorXd _errorWeight;
		
		// Washout position.
		Eigen::VectorXd _washoutPosition;

		// The more the washout factor, the more penalty for the terminal state to be far from the default (washout) position.
		double _washoutFactor;

		Eigen::VectorXd getWashoutState() const;
		void getStateSpaceA(Eigen::MatrixXd& A) const;
		void getStateSpaceB(Eigen::MatrixXd& B) const;
	};
}
