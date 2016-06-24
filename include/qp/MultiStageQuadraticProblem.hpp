#pragma once

#include "MultiStageQuadraticProblemBase.hpp"

#include <Eigen/Dense>

#include <vector>
#include <memory>
#include <ostream>

namespace tmpc
{
	/* Stores data for a multistage QP problem.
	 * Storage format is not explicitly defined and no access to raw data is provided..
	 *
	 *  The problem is stated as following:
	*
	*	min  sum_{ k = 0..nI } z_k'*H_k*z_k + g_k'*z_k
	*	s.t. x_{ k + 1 } = C_k * z_k + c_k				for k = 0..nI - 1
	*            dLow_k <= D_k * z_k <= dUpp_k			for k = 0..nI
	*            zMin_k <= z_k <= zMax_k                for k = 0..nI
	*
	*	where x_k is implicitly defined by z_k = [x_k  u_k] as the first nX variables of z_k
	*
	*	It holds
	*	z_k  \in R^nZ  for k = 0..nI - 1
	*   z_nI \in R*nX
	*
	*	nX < nZ
	*	nU = nZ - nX
	*/
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class MultiStageQuadraticProblem : public MultiStageQuadraticProblemBase<MultiStageQuadraticProblem<NX_, NU_, NC_, NCT_>>
	{
	public:
		typedef unsigned int size_type;

		static unsigned const NX = NX_;
		static unsigned const NU = NU_;
		static unsigned const NC = NC_;
		static unsigned const NCT = NCT_;

		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NU, 1> InputVector;
		typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
		typedef Eigen::Matrix<double, NU, NU> InputInputMatrix;
		typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;
		typedef Eigen::Matrix<double, NC, NX> ConstraintStateMatrix;
		typedef Eigen::Matrix<double, NC, NU> ConstraintInputMatrix;
		typedef Eigen::Matrix<double, NC, 1> ConstraintVector;
		typedef Eigen::Matrix<double, NCT, NX> EndConstraintStateMatrix;
		typedef Eigen::Matrix<double, NCT, 1> EndConstraintVector;

		size_type nT() const { return _nt; }
		static constexpr size_type nX() { return NX; }
		static constexpr size_type nU() { return NU; }
		static constexpr size_type nD() { return NC; }
		static constexpr size_type nDT() { return NCT; }

		StateStateMatrix const& get_Q(size_type i) const { return stage(i, 1).Q; }
		template <typename Matrix> void set_Q(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).Q = val; }

		InputInputMatrix const& get_R(size_type i) const { return stage(i).R; }
		template <typename Matrix> void set_R(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).R = val; }

		StateInputMatrix const& get_S(size_type i) const { return stage(i).S; }
		template <typename Matrix> void set_S(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).S = val; }

		/*
		StateStateMatrix const& get_Q_end() const { return _Q_end; }
		template <typename Matrix> void set_Q_end(Eigen::MatrixBase<Matrix> const& val) { _Q_end = val; }
		*/

		StateVector const& get_q(size_type i) const { return stage(i, 1).q; }
		template <typename Matrix> void set_q(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).q = val; }

		InputVector const& get_r(size_type i) const { return stage(i).r; }
		template <typename Matrix> void set_r(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).r = val; }

		/*
		StateVector const& get_q_end() const { return _q_end;	}
		template <typename Matrix> void set_q_end(Eigen::MatrixBase<Matrix> const& val) { _q_end = val; }
		*/

		StateStateMatrix const& get_A(size_type i) const { return stage(i).A; }
		template <typename Matrix> void set_A(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).A = val; }

		StateInputMatrix const& get_B(size_type i) const { return stage(i).B; }
		template <typename Matrix> void set_B(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).B = val; }

		StateVector const& get_b(size_type i) const { return stage(i).b; }
		template <typename Matrix> void set_b(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).b = val; }

		ConstraintStateMatrix const& get_C(size_type i) const { return stage(i).C; }
		template <typename Matrix> void set_C(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).C = val; }

		ConstraintInputMatrix const& get_D(size_type i) const { return stage(i).D; }
		template <typename Matrix> void set_D(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).D = val; }

		ConstraintStateMatrix const& get_C_end() const { return _C_end; }
		template <typename Matrix> void set_C_end(Eigen::MatrixBase<Matrix> const& val) { _C_end = val; }

		ConstraintVector const& get_d_min(size_type i) const { return stage(i).d_min; }
		template <typename Matrix> void set_d_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).d_min = val; }

		EndConstraintVector const& get_d_end_min() const { return _d_end_min; }
		template <typename Matrix> void set_d_end_min(Eigen::MatrixBase<Matrix> const& val)	{ _d_end_min = val; }

		ConstraintVector const& get_d_max(size_type i) const { return stage(i).d_max; }
		template <typename Matrix> void set_d_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).d_max = val; }

		EndConstraintVector const& get_d_end_max() const { return _d_end_max; }
		template <typename Matrix> void set_d_end_max(Eigen::MatrixBase<Matrix> const& val)	{ _d_end_max = val; }
				
		StateVector const& get_x_min(size_type i) const { return stage(i, 1).x_min; }
		template <typename Matrix> void set_x_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).x_min = val; }

		StateVector const& get_x_max(size_type i) const { return stage(i, 1).x_max; }
		template <typename Matrix> void set_x_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).x_max = val; }

		/*
		StateVector const& get_x_end_min() { return _x_end_min; }
		template <typename Matrix> void set_x_end_min(Eigen::MatrixBase<Matrix> const& val) { _x_end_min = val; }

		StateVector const& get_x_end_max() const { return _x_end_max; }
		template <typename Matrix> void set_x_end_max(Eigen::MatrixBase<Matrix> const& val) { _x_end_max = val; }
		*/
		
		InputVector const& get_u_min(size_type i) const { return stage(i).u_min; }
		template <typename Matrix> void set_u_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).u_min = val; }

		InputVector const& get_u_max(size_type i) const { return stage(i).u_max; }
		template <typename Matrix> void set_u_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).u_max = val; }

		MultiStageQuadraticProblem(size_type nt)
		:	_nt(nt)
		,	_stage(nt + 1)	// only some fields of _stage[nt] are used
		{}

	private:
		struct StageData
		{
			StateStateMatrix Q;
			StateInputMatrix S;
			InputInputMatrix R;
			StateVector q;
			InputVector r;
			StateStateMatrix A;
			StateInputMatrix B;
			StateVector b;
			StateVector x_min;
			StateVector x_max;
			InputVector u_min;
			InputVector u_max;
			ConstraintStateMatrix C;
			ConstraintInputMatrix D;
			ConstraintVector d_min;
			ConstraintVector d_max;
		};

		StageData& stage(size_type i, std::size_t delta = 0)
		{
			assert(i < nT() + delta && i < _stage.size());
			return _stage[i];
		}

		StageData const& stage(size_type i, std::size_t delta = 0) const
		{
			assert(i < nT() + delta && i < _stage.size());
			return _stage[i];
		}

		// Private data members.
		//

		// Number of control intervals.
		std::size_t _nt;

		// Stores stage data
		std::vector<StageData> _stage;

		// 1 matrix of size _Nx x _Nx.
		StateStateMatrix _Q_end;

		// 1 vector of size _Nx
		StateVector _q_end;

		// 1 matrix of size NCT x NX.
		EndConstraintStateMatrix _C_end;

		// 1 vector of size NCT
		EndConstraintVector _d_end_min;

		// 1 vector of size NCT
		EndConstraintVector _d_end_max;

		// 1 vector of size _Nx
		StateVector _x_end_min;

		// 1 vector of size _Nx
		StateVector _x_end_max;
	};
}
