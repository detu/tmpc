#pragma once

namespace tmpc
{
	// Gauss-Newton approximation of a Hessian.
	template <typename CMatrix, typename DMatrix,
		typename QMatrix, typename RMatrix, typename SMatrix>
	void Gauss_Newton_approximation(CMatrix const& C, DMatrix const& D,
			QMatrix& Q, RMatrix& R, SMatrix& S)
	{
		// H = G^T G
		//   = [Q S
		//      S R]
		//

		Q = transpose(C) * C;
		R = transpose(D) * D;
		S = transpose(C) * D;
	}

	// Gauss-Newton approximation of a Hessian. The gradient is also computed.
	template <typename ResidualVector, typename CMatrix, typename DMatrix,
		typename QMatrix, typename RMatrix, typename SMatrix, typename StateGradientVector, typename InputGradientVector>
	void Gauss_Newton_approximation(ResidualVector const& res, CMatrix const& C, DMatrix const& D,
			QMatrix& Q, RMatrix& R, SMatrix& S, StateGradientVector& q, InputGradientVector& r)
	{
		// H = G^T G
		//   = [Q S
		//      S R]
		//

		Q = transpose(C) * C;
		R = transpose(D) * D;
		S = transpose(C) * D;

		// g = 2 * (y_bar - y_hat)^T * W * G
		// g = [q; r]
		q = transpose(res) * C;
		r = transpose(res) * D;
	}

	// Gauss-Newton approximation of a Hessian with weighting matrix.
	template <typename ResidualVector, typename WeightingMatrix, typename CMatrix, typename DMatrix,
		typename QMatrix, typename RMatrix, typename SMatrix, typename StateGradientVector, typename InputGradientVector>
	void Gauss_Newton_approximation(ResidualVector const& res, WeightingMatrix const& W, CMatrix const& C, DMatrix const& D,
			QMatrix& Q, RMatrix& R, SMatrix& S, StateGradientVector& q, InputGradientVector& r)
	{
		// H = G^T W G
		//   = [Q S
		//      S R]
		//

		Q = transpose(C) * W * C;
		R = transpose(D) * W * D;
		S = transpose(C) * W * D;

		// g = 2 * (y_bar - y_hat)^T * W * G
		// g = [q; r]
		q = transpose(res) * W * C;
		r = transpose(res) * W * D;
	}
}
