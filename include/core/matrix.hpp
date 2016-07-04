/*
 * Common matrix types and operations.
 */

#pragma once

#include <Eigen/Dense>

#include <limits>

namespace tmpc
{
	template <typename Matrix>
	typename Matrix::ConstantReturnType signaling_nan(typename Matrix::Index M, typename Matrix::Index N)
	{
		return Matrix::Constant(M, N, std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	typename Matrix::ConstantReturnType signaling_nan(typename Matrix::Index N)
	{
		return Matrix::Constant(N, std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	typename Matrix::ConstantReturnType signaling_nan()
	{
		return Matrix::Constant(std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	typename Matrix::ConstTransposeReturnType transpose(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.transpose();
	}

	template <std::size_t N, typename Matrix>
	decltype(auto) middle_rows(Eigen::MatrixBase<Matrix> const& m, std::size_t first_row)
	{
		return m.template middleRows<N>(first_row);
	}
}
