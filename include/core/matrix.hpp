/*
 * Common matrix types and operations.
 */

#pragma once

#include <Eigen/Dense>

#include <limits>
#include <type_traits>

namespace tmpc
{
	template <typename Matrix>
	typename Matrix::ConstantReturnType constant(typename Matrix::Index M, typename Matrix::Index N, typename Matrix::Scalar val)
	{
		return Matrix::Constant(M, N, val);
	}

	template <typename Matrix>
	typename Matrix::ConstantReturnType constant(typename Matrix::Index N, typename Matrix::Scalar val)
	{
		return Matrix::Constant(N, val);
	}

	template <typename Matrix>
	typename Matrix::ConstantReturnType constant(typename Matrix::Scalar val)
	{
		return Matrix::Constant(val);
	}

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

	template <std::size_t N, typename Matrix>
	decltype(auto) middle_rows(Eigen::MatrixBase<Matrix>& m, std::size_t first_row)
	{
		return m.template middleRows<N>(first_row);
	}

	template<unsigned N, typename Matrix>
	decltype(auto) top_rows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template topRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) top_rows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) bottom_rows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) bottom_rows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) left_cols(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template leftCols<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) left_cols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template leftCols<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) right_cols(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template rightCols<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) right_cols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template rightCols<N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	decltype(auto) top_left_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topLeftCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	decltype(auto) top_right_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topRightCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	decltype(auto) bottom_left_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomLeftCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	decltype(auto) bottom_right_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomRightCorner<M, N>();
	}

	template <typename Matrix>
	std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::IdentityReturnType> identity()
	{
		return Matrix::Identity();
	}

	template <typename Matrix>
	std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::Index> constexpr rows()
	{
		return Matrix::RowsAtCompileTime;
	}

	template <typename Matrix>
	std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::Index> constexpr cols()
	{
		return Matrix::ColsAtCompileTime;
	}

	template <typename Matrix>
	typename Matrix::EvalReturnType eval(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.eval();
	}

	template <typename Matrix>
	// Eigen::NoAlias< Matrix, Eigen::MatrixBase<Matrix> > does not work as a return type
	// for the reason I which don't understand:
	// error: template argument for template template parameter must be a class template or type alias template
    // Eigen::NoAlias< Matrix, Eigen::MatrixBase<Matrix> > noalias(Eigen::MatrixBase<Matrix> const& m)
    //                         ^
	// Eigen::NoAlias< Matrix, Eigen::MatrixBase<Matrix> > noalias(Eigen::MatrixBase<Matrix> const& m)
	decltype(auto) noalias(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.noalias();
	}

	template <typename Matrix>
	decltype(auto) as_diagonal(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.asDiagonal();
	}

	template <typename Matrix>
	decltype(auto) diagonal(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.diagonal();
	}

	template <typename Matrix>
	decltype(auto) diagonal(Eigen::MatrixBase<Matrix>& m)
	{
		return m.diagonal();
	}
}
