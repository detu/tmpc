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

	static auto constexpr NX = NX_;
	static auto constexpr NU = NU_;
	static auto constexpr NY = NY_;
	static auto constexpr NW = NW_;
	static auto constexpr NP = NP_;
	static auto constexpr NC = NC_;
	static auto constexpr NCT = NCT_;

	template <unsigned M, unsigned N>
	using Matrix = Eigen::Matrix<Scalar, M, N, EigenOptions>;

	typedef Matrix<NP, 1> ParameterVector;
	typedef Matrix<NU, 1> InputVector;
	typedef Matrix<NW, 1> DisturbanceVector;
	typedef Matrix<NX, 1> StateVector;
	typedef Matrix<NY, 1> OutputVector;
	typedef Matrix<NC, 1> ConstraintVector;
	typedef Matrix<NX, NX> StateStateMatrix;
	typedef Matrix<NX, NU> StateInputMatrix;
	typedef Matrix<NU, NU> InputInputMatrix;
	typedef Matrix<NX, NW> StateDisturbanceMatrix;
	typedef Matrix<NC, NU> ConstraintInputMatrix;
	typedef Matrix<NY, NX> OutputStateMatrix;
	typedef Matrix<NY, NU> OutputInputMatrix;
	typedef Matrix<NY, NY> OutputOutputMatrix;
	typedef Matrix<NW, NW> DisturbanceDisturbanceMatrix;

	/*
	template <typename Matrix>
	using LLT = Eigen::LLT<Matrix>;
	*/
};

}
