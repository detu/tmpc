#pragma once

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
template <typename K, typename D>
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

	typedef typename K::template Vector<NP> ParameterVector;
	typedef typename K::template Vector<NU> InputVector;
	typedef typename K::template Vector<NW> DisturbanceVector;
	typedef typename K::template Vector<NX> StateVector;
	typedef typename K::template Vector<NY> OutputVector;
	typedef typename K::template Vector<NC> ConstraintVector;
	typedef typename K::template Vector<NCT> TerminalConstraintVector;

	typedef typename K::template Matrix<NX, NX> StateStateMatrix;
	typedef typename K::template Matrix<NX, NU> StateInputMatrix;
	typedef typename K::template Matrix<NX, NW> StateDisturbanceMatrix;

	typedef typename K::template Matrix<NU, NU> InputInputMatrix;

	typedef typename K::template Matrix<NW, NW> DisturbanceDisturbanceMatrix;

	typedef typename K::template Matrix<NC, NU> ConstraintInputMatrix;
	typedef typename K::template Matrix<NC, NX> ConstraintStateMatrix;
	typedef typename K::template Matrix<NC, NC> ConstraintConstraintMatrix;

	typedef typename K::template Matrix<NCT, NU> TerminalConstraintInputMatrix;
	typedef typename K::template Matrix<NCT, NX> TerminalConstraintStateMatrix;

	typedef typename K::template Matrix<NY, NX> OutputStateMatrix;
	typedef typename K::template Matrix<NY, NU> OutputInputMatrix;
	typedef typename K::template Matrix<NY, NW> OutputDisturbanceMatrix;
	typedef typename K::template Matrix<NY, NY> OutputOutputMatrix;
};

}	// namespace tmpc
