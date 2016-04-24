#pragma once

#include <ModelPredictiveController.hpp>
#include <CyberMotionOCP.hpp>

namespace mpmc
{
	class MotionCueingController : public camels::ModelPredictiveController<CyberMotionOCP>
	{
	public:
		typedef CyberMotionOCP::StateInputVector StateInputVector;
		typedef CyberMotionOCP::OutputVector OutputVector;
		typedef CyberMotionOCP::OutputJacobianMatrix OutputJacobianMatrix;

		MotionCueingController(CyberMotionOCP const& ocp, double sample_time);

	private:
		void PathConstraints(unsigned i, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const override;
		void TerminalConstraints(const Eigen::VectorXd& x, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const override;

		// Private data members.
		const unsigned _Nq;
		const unsigned _Ny;

		template<class Vector1, class Matrix, class Vector2>
		void SRConstraints(const Eigen::MatrixBase<Vector1>& x, Eigen::MatrixBase<Matrix>& D, Eigen::MatrixBase<Vector2>& d_min, Eigen::MatrixBase<Vector2>& d_max) const;
	};
}
