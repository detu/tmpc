#pragma once

#include <Eigen/Dense>

namespace tmpc
{

/**
 * \brief Class defining problem dimensions, vector and matrix types,
 * and algebraic operations using Eigen3 library.
 *
 * \tparam Scalar_ scalar type used for matrices and vectors
 * \tparam NX_ number of states
 * \tparam NU_ number of inputs
 * \tparam NW_ number of unmeasured disturbances
 * \tparam NY_ number of measured outputs
 * \tparam NP_ number of system parameters
 * \tparam NC_ number of path constraints
 * \tparam NCT_ number of terminal path constraints
 * \tparam EigenOptions options for Eigen controlling matrix layout (ColMajor/RowMajor) and alignment.
 */
template <typename Scalar_, unsigned NX_, unsigned NU_, unsigned NW_,
	unsigned NY_, unsigned NP_, unsigned NC_, unsigned NCT_, int EigenOptions = Eigen::ColMajor>
class EigenKernel
{
public:
	EigenKernel() = delete;

	typedef Scalar_ Scalar;
	typedef Eigen::Index size_t;

	static size_t constexpr NX = NX_;
	static size_t constexpr NU = NU_;
	static size_t constexpr NY = NY_;
	static size_t constexpr NW = NW_;
	static size_t constexpr NP = NP_;
	static size_t constexpr NC = NC_;
	static size_t constexpr NCT = NCT_;

	template <unsigned M, unsigned N>
	using Matrix = Eigen::Matrix<Scalar, M, N, EigenOptions>;

	template <unsigned M>
	using Vector = Eigen::Matrix<Scalar, M, 1, EigenOptions>;

	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, EigenOptions> DynamicMatrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, EigenOptions> DynamicVector;

	typedef Vector<NP> ParameterVector;
	typedef Vector<NU> InputVector;
	typedef Vector<NW> DisturbanceVector;
	typedef Vector<NX> StateVector;
	typedef Vector<NY> OutputVector;
	typedef Vector<NC> ConstraintVector;
	typedef Matrix<NX, NX> StateStateMatrix;
	typedef Matrix<NX, NU> StateInputMatrix;
	typedef Matrix<NU, NU> InputInputMatrix;
	typedef Matrix<NX, NW> StateDisturbanceMatrix;
	typedef Matrix<NC, NU> ConstraintInputMatrix;
	typedef Matrix<NY, NX> OutputStateMatrix;
	typedef Matrix<NY, NU> OutputInputMatrix;
	typedef Matrix<NY, NW> OutputDisturbanceMatrix;
	typedef Matrix<NY, NY> OutputOutputMatrix;
	typedef Matrix<NW, NW> DisturbanceDisturbanceMatrix;

	/*
	template <typename Matrix>
	using LLT = Eigen::LLT<Matrix>;
	*/

	static typename DynamicMatrix::ConstantReturnType constant(size_t M, size_t N, Scalar val)
	{
		return DynamicMatrix::Constant(M, N, val);
	}

	static typename DynamicVector::ConstantReturnType constant(size_t N, Scalar val)
	{
		return DynamicVector::Constant(N, val);
	}

	template <typename Matrix>
	static typename Matrix::ConstantReturnType constant(Scalar val)
	{
		return Matrix::Constant(val);
	}

	static typename DynamicMatrix::ConstantReturnType zero(size_t M, size_t N)
	{
		return DynamicMatrix::Zero(M, N);
	}

	static typename DynamicVector::ConstantReturnType zero(size_t N)
	{
		return DynamicVector::Zero(N);
	}

	template <typename Matrix>
	static typename Matrix::ConstantReturnType zero()
	{
		return Matrix::Zero();
	}

	template <typename Matrix>
	static void set_zero(Eigen::MatrixBase<Matrix>& m)
	{
		m.setZero();
	}

	template <typename Matrix>
	static typename Matrix::ConstantReturnType signaling_nan(typename Matrix::Index M, typename Matrix::Index N)
	{
		return Matrix::Constant(M, N, std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	static typename Matrix::ConstantReturnType signaling_nan(typename Matrix::Index N)
	{
		return Matrix::Constant(N, std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	static typename Matrix::ConstantReturnType signaling_nan()
	{
		return Matrix::Constant(std::numeric_limits<typename Matrix::Scalar>::signaling_NaN());
	}

	template <typename Matrix>
	static typename Matrix::ConstTransposeReturnType transpose(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.transpose();
	}

	template <typename Matrix>
	static decltype(auto) row(Eigen::MatrixBase<Matrix> const& m, typename Matrix::Index i)
	{
		return m.row(i);
	}

	template <typename Matrix>
	static decltype(auto) row(Eigen::MatrixBase<Matrix>& m, typename Matrix::Index i)
	{
		return m.row(i);
	}

	template <typename Matrix>
	static decltype(auto) col(Eigen::MatrixBase<Matrix> const& m, typename Matrix::Index i)
	{
		return m.col(i);
	}

	template <typename Matrix>
	static decltype(auto) col(Eigen::MatrixBase<Matrix>& m, typename Matrix::Index i)
	{
		return m.col(i);
	}

	template <std::size_t N, typename Matrix>
	static decltype(auto) middle_rows(Eigen::MatrixBase<Matrix> const& m, std::size_t first_row)
	{
		return m.template middleRows<N>(first_row);
	}

	template <std::size_t N, typename Matrix>
	static decltype(auto) middle_rows(Eigen::MatrixBase<Matrix>& m, std::size_t first_row)
	{
		return m.template middleRows<N>(first_row);
	}

	template <unsigned N, typename Matrix>
	static decltype(auto) top_rows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template topRows<N>();
	}

	template <unsigned N, typename Matrix>
	static decltype(auto) top_rows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topRows<N>();
	}

	template <typename Matrix>
	static decltype(auto) top_rows(Eigen::MatrixBase<Matrix>& m, size_t N)
	{
		return m.template topRows(N);
	}

	template <typename Matrix>
	static decltype(auto) top_rows(Eigen::MatrixBase<Matrix> const& m, size_t N)
	{
		return m.template topRows(N);
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) bottom_rows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) bottom_rows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) left_cols(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template leftCols<N>();
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) left_cols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template leftCols<N>();
	}

	template<typename Matrix>
	static decltype(auto) left_cols(Eigen::MatrixBase<Matrix>& m, size_t N)
	{
		return m.template leftCols(N);
	}

	template<typename Matrix>
	static decltype(auto) left_cols(Eigen::MatrixBase<Matrix> const& m, size_t N)
	{
		return m.template leftCols(N);
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) right_cols(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template rightCols<N>();
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) right_cols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template rightCols<N>();
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) middle_cols(Eigen::MatrixBase<Matrix>& m, size_t j)
	{
		return m.template middleCols<N>(j);
	}

	template<unsigned N, typename Matrix>
	static decltype(auto) middle_cols(Eigen::MatrixBase<Matrix> const& m, size_t j)
	{
		return m.template middleCols<N>(j);
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) top_left_corner(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template topLeftCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) top_left_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topLeftCorner<M, N>();
	}

	template<typename Matrix>
	static decltype(auto) top_left_corner(Eigen::MatrixBase<Matrix>& m, size_t M, size_t N)
	{
		return m.template topLeftCorner(M, N);
	}

	template<typename Matrix>
	static decltype(auto) top_left_corner(Eigen::MatrixBase<Matrix> const& m, size_t M, size_t N)
	{
		return m.template topLeftCorner(M, N);
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) top_right_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topRightCorner<M, N>();
	}

	template <typename Matrix>
	static decltype(auto) top_right_corner(Eigen::MatrixBase<Matrix>& m, size_t M, size_t N)
	{
		return m.template topRightCorner(M, N);
	}

	template <typename Matrix>
	static decltype(auto) top_right_corner(Eigen::MatrixBase<Matrix> const& m, size_t M, size_t N)
	{
		return m.template topRightCorner(M, N);
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) bottom_left_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomLeftCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) bottom_right_corner(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template bottomRightCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) bottom_right_corner(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomRightCorner<M, N>();
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) block(Eigen::MatrixBase<Matrix>& m, size_t i, size_t j)
	{
		return m.template block<M, N>(i, j);
	}

	template<unsigned M, unsigned N, typename Matrix>
	static decltype(auto) block(Eigen::MatrixBase<Matrix> const& m, size_t i, size_t j)
	{
		return m.template block<M, N>(i, j);
	}

	template <typename Matrix>
	static std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::IdentityReturnType> identity()
	{
		return Matrix::Identity();
	}

	static typename DynamicMatrix::IdentityReturnType identity(size_t M, size_t N)
	{
		return DynamicMatrix::Identity(M, N);
	}

	template <typename Matrix>
	static std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::Index> constexpr rows()
	{
		return Matrix::RowsAtCompileTime;
	}

	template <typename Matrix>
	static size_t rows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.rows();
	}

	template <typename Matrix>
	static std::enable_if_t<std::is_base_of<Eigen::MatrixBase<Matrix>, Matrix>::value, typename Matrix::Index> constexpr cols()
	{
		return Matrix::ColsAtCompileTime;
	}

	template <typename Matrix>
	static size_t cols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.cols();
	}

	template <typename Matrix>
	static typename Matrix::EvalReturnType eval(Eigen::MatrixBase<Matrix> const& m)
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
	static decltype(auto) noalias(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.noalias();
	}

	template <typename Matrix>
	static decltype(auto) as_diagonal(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.asDiagonal();
	}

	template <typename Matrix>
	static decltype(auto) diagonal(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.diagonal();
	}

	template <typename Matrix>
	static decltype(auto) diagonal(Eigen::MatrixBase<Matrix>& m)
	{
		return m.diagonal();
	}

	template <typename Matrix>
	static decltype(auto) squared_norm(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.squaredNorm();
	}

	template <typename Matrix>
	static decltype(auto) norm_2(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.norm();
	}

	template <typename Matrix>
	static decltype(auto) inverse(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.inverse();
	}

	template <typename Matrix>
	static decltype(auto) selfadjoint_view_upper(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template selfadjointView<Eigen::Upper>();
	}

	template <typename Matrix>
	static decltype(auto) selfadjoint_view_upper(Eigen::MatrixBase<Matrix>&& m)
	{
		return m.template selfadjointView<Eigen::Upper>();
	}

	template <typename Matrix>
	static decltype(auto) selfadjoint_view_upper(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template selfadjointView<Eigen::Upper>();
	}

	template <typename Matrix>
	static decltype(auto) triangular_view_upper(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template triangularView<Eigen::Upper>();
	}

	template <typename Matrix>
	static decltype(auto) triangular_view_upper(Eigen::MatrixBase<Matrix>&& m)
	{
		return m.template triangularView<Eigen::Upper>();
	}

	template <typename Matrix>
	static decltype(auto) triangular_view_upper(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template triangularView<Eigen::Upper>();
	}
};

}
