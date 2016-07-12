#pragma once

namespace tmpc
{
	// Gauss-Newton approximation of a Hessian.
	template <typename ResidualVector, typename WeightingMatrix, typename CMatrix, typename DMatrix,
		typename QMatrix, typename RMatrix, typename SMatrix, typename StateGradientVector, typename InputGradientVector>
	void Gauss_Newton_approximation(ResidualVector const& res, WeightingMatrix const& W, CMatrix const& C, DMatrix const& D,
			QMatrix& Q, RMatrix& R, SMatrix& S, StateGradientVector& q, InputGradientVector& r)
	{
		// H = G^T W G + \mu I
		//   = [Q S
		//      S R]
		//

		/*
		 *  For the reason that I don't understand, Eigen does not want to digest this:
		 *  Q = tmpc::transpose(C) * W * C + tmpc::as_diagonal(_stateErrorWeight);
		 *
		 *  It says:
			/home/kotlyar/Documents/MPI/deathstar/src/deathstar_mpc/CableRobotOCP.hpp:72:35: error: invalid operands to binary expression ('const Product<Eigen::Product<Eigen::Transpose<const Eigen::Block<Eigen::Matrix<double, 9, 21, 0, 9, 21>, 9, 13, true> >, Eigen::DiagonalWrapper<const Eigen::CwiseUnaryOp<Eigen::internal::scalar_abs2_op<double>, const Eigen::Matrix<double, 9, 1, 0, 9, 1> > >, 1>, Eigen::Block<Eigen::Matrix<double, 9, 21, 0, 9, 21>, 9, 13, true> >' and 'const Eigen::DiagonalWrapper<const Eigen::Matrix<double, 13, 1, 0, 13, 1> >')
					Q = tmpc::transpose(C) * W * C + tmpc::as_diagonal(_stateErrorWeight);
						~~~~~~~~~~~~~~~~~~~~~~~~~~ ^ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			/usr/local/include/eigen3/Eigen/src/Core/arch/CUDA/Half.h:189:38: note: candidate function not viable: no known conversion from 'const Product<Eigen::Product<Eigen::Transpose<const Eigen::Block<Eigen::Matrix<double, 9, 21, 0, 9, 21>, 9, 13, true> >, Eigen::DiagonalWrapper<const Eigen::CwiseUnaryOp<Eigen::internal::scalar_abs2_op<double>, const Eigen::Matrix<double, 9, 1, 0, 9, 1> > >, 1>, Eigen::Block<Eigen::Matrix<double, 9, 21, 0, 9, 21>, 9, 13, true> >' to 'const Eigen::half' for 1st argument
			static inline EIGEN_DEVICE_FUNC half operator + (const half& a, const half& b) {
												 ^
			/usr/local/include/eigen3/Eigen/src/Core/../plugins/CommonCwiseBinaryOps.h:27:28: note: candidate template ignored: could not match 'MatrixBase' against 'DiagonalWrapper'
			EIGEN_MAKE_CWISE_BINARY_OP(operator+,internal::scalar_sum_op)

			Therefore, I need to write it in two lines using +=.
		 */
		Q = tmpc::transpose(C) * W * C;
		R = tmpc::transpose(D) * W * D;
		S = tmpc::transpose(C) * W * D;

		// g = 2 * (y_bar - y_hat)^T * W * G
		// g = [q; r]
		q = tmpc::transpose(res) * W * C;
		r = tmpc::transpose(res) * W * D;
	}
}
