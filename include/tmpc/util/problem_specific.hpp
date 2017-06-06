#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc {

/**
 * \brief Defines problem-specific matrix and vector types.
 *
 * \tparam <K> Class defining matrix and vector types.
 * \tparam <D> Class defining problem sizes. D must define the following constants:
 *   - NX number of states
 *   - NU number of inputs
 *   - NW number of unmeasured disturbances
 *   - NY number of measured outputs
 *   - NP number of system parameters
 *   - NC number of path constraints
 *   - NCT number of terminal path constraints
 *
 *   Hint: can be conviniently used as a mixin-type.
 */
template <typename Scalar, typename D>
class ProblemSpecific
{
public:
	static auto constexpr NX = D::NX;
	static auto constexpr NU = D::NU;
	static auto constexpr NY = D::NY;
	static auto constexpr NW = D::NW;
	static auto constexpr NP = D::NP;
	static auto constexpr NC = D::NC;
	static auto constexpr NCT = D::NCT;

	typedef StaticVector<Scalar, NP> ParameterVector;
	typedef StaticVector<Scalar, NU> InputVector;
	typedef StaticVector<Scalar, NW> DisturbanceVector;
	typedef StaticVector<Scalar, NX> StateVector;
	typedef StaticVector<Scalar, NY> OutputVector;
	typedef StaticVector<Scalar, NC> ConstraintVector;
	typedef StaticVector<Scalar, NCT> TerminalConstraintVector;

	typedef StaticMatrix<Scalar, NX, NX> StateStateMatrix;
	typedef StaticMatrix<Scalar, NX, NU> StateInputMatrix;
	typedef StaticMatrix<Scalar, NX, NW> StateDisturbanceMatrix;

	typedef StaticMatrix<Scalar, NU, NU> InputInputMatrix;

	typedef StaticMatrix<Scalar, NW, NW> DisturbanceDisturbanceMatrix;

	typedef StaticMatrix<Scalar, NC, NU> ConstraintInputMatrix;
	typedef StaticMatrix<Scalar, NC, NX> ConstraintStateMatrix;
	typedef StaticMatrix<Scalar, NC, NC> ConstraintConstraintMatrix;

	typedef StaticMatrix<Scalar, NCT, NU> TerminalConstraintInputMatrix;
	typedef StaticMatrix<Scalar, NCT, NX> TerminalConstraintStateMatrix;

	typedef StaticMatrix<Scalar, NY, NX> OutputStateMatrix;
	typedef StaticMatrix<Scalar, NY, NU> OutputInputMatrix;
	typedef StaticMatrix<Scalar, NY, NW> OutputDisturbanceMatrix;
	typedef StaticMatrix<Scalar, NY, NY> OutputOutputMatrix;
};

}	// namespace tmpc
