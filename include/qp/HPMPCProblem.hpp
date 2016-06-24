#pragma once

#include "MultiStageQuadraticProblemBase.hpp"

#include <Eigen/Dense>

namespace tmpc
{
	template< unsigned M, unsigned N >
	struct HPMPCMatrixType
	{
		typedef Eigen::Matrix<double, M, N, Eigen::RowMajor> type;
	};

	template< unsigned M >
	struct HPMPCMatrixType<M, 1>
	{
		typedef Eigen::Matrix<double, M, 1> type;
	};

	template< unsigned N >
	struct HPMPCMatrixType<1, N>
	{
		typedef Eigen::Matrix<double, 1, N, Eigen::RowMajor> type;
	};

	template<>
	struct HPMPCMatrixType<1, 1>
	{
		typedef Eigen::Matrix<double, 1, 1> type;
	};

	template< unsigned NX, unsigned NU >
	using HPMPCMatrix = typename HPMPCMatrixType<NX, NU>::type;

	/* Stores data for a multistage QP problem.
	 * The storage format is what is expected by the c_order_d_ip_ocp_hard_tv() function from HPMPC.
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
	class HPMPCProblem : public MultiStageQuadraticProblemBase<HPMPCProblem<NX_, NU_, NC_, NCT_>>
	{
	public:
		typedef unsigned int size_type;

		static unsigned const NX = NX_;
		static unsigned const NU = NU_;
		static unsigned const NZ = NX + NU;
		static unsigned const NC = NC_;
		static unsigned const NCT = NCT_;

		typedef HPMPCMatrix<NX ,  1> StateVector;
		typedef HPMPCMatrix<NU ,  1> InputVector;
		typedef HPMPCMatrix<NZ ,  1> StateInputVector;
		typedef HPMPCMatrix<NZ , NZ> StageHessianMatrix;
		typedef HPMPCMatrix<NX , NX> EndStageHessianMatrix;
		typedef HPMPCMatrix<NZ ,  1> StageGradientVector;
		typedef HPMPCMatrix<NX ,  1> EndStageGradientVector;
		typedef HPMPCMatrix<NX , NZ> InterStageMatrix;
		typedef HPMPCMatrix<NC , NZ> StageConstraintMatrix;
		typedef HPMPCMatrix<NC ,  1> StageConstraintVector;
		typedef HPMPCMatrix<NCT, NX> EndStageConstraintMatrix;
		typedef HPMPCMatrix<NCT,  1> EndStageConstraintVector;

		typedef HPMPCMatrix<NU, NU> HPMPC_RMatrix;
		typedef HPMPCMatrix<NU, NX> HPMPC_SMatrix;
		typedef HPMPCMatrix<NX, NX> HPMPC_QMatrix;
		typedef HPMPCMatrix<NX, NX> HPMPC_AMatrix;
		typedef HPMPCMatrix<NX, NU> HPMPC_BMatrix;
		typedef HPMPCMatrix<NC, NX> HPMPC_CMatrix;
		typedef HPMPCMatrix<NC, NU> HPMPC_DMatrix;

		size_type nT() const { return _nt; }
		static constexpr size_type nX() { return NX; }
		static constexpr size_type nZ() { return NZ; }
		static constexpr size_type nU() { return NU; }
		static constexpr size_type nD() { return NC; }
		static constexpr size_type nDT() { return NCT; }

		template<class Matrix> void set_Q(std::size_t i, Eigen::MatrixBase<Matrix> const& Q) { stage(i, 1)._Q = Q; }
		HPMPC_QMatrix const& get_Q(std::size_t i) const { return stage(i, 1)._Q; }

		template<class Matrix> void set_R(std::size_t i, Eigen::MatrixBase<Matrix> const& R) { stage(i)._R = R; }
		HPMPC_RMatrix const& get_R(std::size_t i) const { return stage(i)._R; }

		// HPMPC "S" size is NUxNX, whereas the MultiStageQuadraticProblem interface assumes NXxNU,
		// hence the transposes.
		template<class Matrix> void set_S(std::size_t i, Eigen::MatrixBase<Matrix> const& S) { stage(i)._S = S.transpose(); }
		decltype(auto) get_S(std::size_t i) const { return stage(i)._S.transpose(); }

		template<class Matrix> void set_q(std::size_t i, Eigen::MatrixBase<Matrix> const& q) { stage(i, 1)._q = q; }
		StateVector const& get_q(std::size_t i) const { return stage(i, 1)._q; }

		template<class Matrix> void set_r(std::size_t i, Eigen::MatrixBase<Matrix> const& r) { stage(i)._r = r; }
		InputVector const& get_r(std::size_t i) const { return stage(i)._r; }

		template<class Matrix> void set_A(std::size_t i, Eigen::MatrixBase<Matrix> const& A) { stage(i)._A = A; }
		HPMPC_AMatrix const& get_A(std::size_t i) const { return stage(i)._A; }

		template<class Matrix> void set_B(std::size_t i, Eigen::MatrixBase<Matrix> const& B) { stage(i)._B = B; }
		HPMPC_BMatrix const& get_B(std::size_t i) const { return stage(i)._B; }

		template<class Matrix> void set_b(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._b = val; }
		StateVector const& get_b(std::size_t i) const { return stage(i)._b; }

		template <typename Matrix> void set_C(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._C = val; }
		HPMPC_CMatrix const& get_C(std::size_t i) const { return stage(i)._C; }

		template <typename Matrix> void set_D(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._D = val; }
		HPMPC_DMatrix const& get_D(std::size_t i) const { return stage(i)._D; }

		template <typename Matrix> void set_d_min(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._dMin = val; }
		StageConstraintVector const& get_d_min(size_type i) const { return stage(i)._dMin; }

		template <typename Matrix> void set_d_max(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._dMax = val; }
		StageConstraintVector const& get_d_max(size_type i) const { return stage(i)._dMax; }

		template <typename Matrix> void set_C_end(Eigen::MatrixBase<Matrix> const& val) { _C_end = val; }
		EndStageConstraintMatrix const& get_C_end() const { return _C_end; }

		template <typename Matrix> void set_d_end_min(Eigen::MatrixBase<Matrix> const& val) { _d_end_min = val; }
		EndStageConstraintVector const& get_d_end_min() const	{ return _d_end_min; }

		template <typename Matrix> void set_d_end_max(Eigen::MatrixBase<Matrix> const& val) { _d_end_max = val; }
		EndStageConstraintVector const& get_d_end_max() const	{ return _d_end_max; }

		template<class Matrix> void set_x_min(std::size_t i, Eigen::MatrixBase<Matrix> const& val)	{
			stage(i, 1)._lb.template middleRows<NX>(i < nT() ? NU : 0) = val;
		}

		decltype(auto) get_x_min(std::size_t i) const {
			return stage(i, 1)._lb.template middleRows<NX>(i < nT() ? NU : 0);
		}

		template<class Matrix> void set_x_max(std::size_t i, Eigen::MatrixBase<Matrix> const& val)	{
			stage(i, 1)._ub.template middleRows<NX>(i < nT() ? NU : 0) = val;
		}

		decltype(auto) get_x_max(std::size_t i) const {
			return stage(i, 1)._ub.template middleRows<NX>(i < nT() ? NU : 0);
		}

		template<class Matrix> void set_u_min(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._lb.template topRows<NU>() = val; }
		decltype(auto) get_u_min(std::size_t i) const { return stage(i)._lb.template topRows<NU>(); }

		template<class Matrix> void set_u_max(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._ub.template topRows<NU>() = val; }
		decltype(auto) get_u_max(std::size_t i) const { return stage(i)._ub.template topRows<NU>(); }

		bool is_x0_equality_constrained() const { return get_x_min(0) == get_x_max(0); }

		// ******************************************************
		//                HPMPC raw data interface.
		//
		// The prefixes before _data() correspond to the names of
		// the argument to c_order_d_ip_ocp_hard_tv().
		// ******************************************************
		double const * const * A_data () const { return _A .data(); }
		double const * const * B_data () const { return _B .data(); }

		double const * const * b_data () const
		{
			if (is_x0_equality_constrained())
			{
				_b0 = _stage[0]._A * get_x_min(0) + _stage[0]._b;
				_b[0] = _b0.data();
			}
			else
			{
				_b[0] = _stage[0]._b.data();
			}

			return _b.data();
		}

		double const * const * Q_data () const { return _Q .data(); }
		double const * const * S_data () const { return _S .data(); }
		double const * const * R_data () const { return _R .data(); }
		double const * const * q_data () const { return _q .data(); }

		double const * const * r_data () const
		{
			if (is_x0_equality_constrained())
			{
				_r0 = _stage[0]._S * get_x_min(0) + _stage[0]._r;
				_r[0] = _r0.data();
			}
			else
			{
				_r[0] = _stage[0]._r.data();
			}

			return _r .data();
		}

		double const * const * lb_data() const { return _lb.data(); }
		double const * const * ub_data() const { return _ub.data(); }
		double const * const * C_data () const { return _C .data(); }
		double const * const * D_data () const { return _D .data(); }
		double const * const * lg_data() const { return _lg.data(); }	// TODO: beware of x[0] being equality constrained!
		double const * const * ug_data() const { return _ug.data(); }	// TODO: beware of x[0] being equality constrained!

		HPMPCProblem(size_type nt)
		:	_nt(nt)
		,	_stage(nt + 1)	// not all fields of _stage[nt] are used
		,	_A (nt)
		,	_B (nt)
		,	_b (nt)
		,	_Q (nt + 1)
		,	_S (nt + 1)
		,	_R (nt + 1)
		,	_q (nt + 1)
		,	_r (nt + 1)
		,	_lb(nt + 1)
		,	_ub(nt + 1)
		,	_C (nt + 1)
		,	_D (nt + 1)
		,	_lg(nt + 1)
		,	_ug(nt + 1)
		{
			// Filling out pointer arrays for HPMPC
			for (std::size_t i = 0; i <= nt; ++i)
			{
				if (i < nt)
				{
					_A[i] = _stage[i]._A.data();
					_B[i] = _stage[i]._B.data();
					_b[i] = _stage[i]._b.data();
				}

				_Q [i] = _stage[i]._Q .data();
				_S [i] = i < nt ? _stage[i]._S .data() : nullptr;
				_R [i] = i < nt ? _stage[i]._R .data() : nullptr;
				_q [i] = _stage[i]._q .data();
				_r [i] = i < nt ? _stage[i]._r .data() : nullptr;
				_lb[i] = _stage[i]._lb.data();
				_ub[i] = _stage[i]._ub.data();
				_C [i] = i < nt ? _stage[i]._C .data() : _C_end.data();
				_D [i] = i < nt ? _stage[i]._D .data() : nullptr;
				_lg[i] = i < nt ? _stage[i]._dMin.data() : _d_end_min.data();
				_ug[i] = i < nt ? _stage[i]._dMax.data() : _d_end_max.data();
			}
		}

		// *************************************************************************
		// OBSOLETE INTERFACE
		// *************************************************************************
		/*
		template<class Matrix> friend void set_Hend(HPMPCProblem& p, Eigen::MatrixBase<Matrix> const& Hend) { p._Hend = Hend; }
		friend EndStageHessianMatrix const& get_Hend(HPMPCProblem const& p) { return p._Hend; }

		template<class Matrix> friend void set_gend(HPMPCProblem& p, Eigen::MatrixBase<Matrix>& val) { p._gend = val; }
		friend EndStageGradientVector const& get_gend(HPMPCProblem const& p) { return p._gend;	}

		EndStageConstraintMatrix const& get_Dend() const { return _C_end; }
		StageConstraintVector const& get_dMin(size_type i) const { return stage(i)._dMin; }
		EndStageConstraintVector const& get_dendMin() const	{ return _d_end_min; }
		StageConstraintVector const& get_dMax(size_type i) const { return stage(i)._dMax; }
		EndStageConstraintVector const& get_dendMax() const	{ return _d_end_max; }

		template<class Matrix> void set_xMin(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._lb.template bottomRows<NX>() = val; }
		decltype(auto) get_xMin(std::size_t i) const { return stage(i)._lb.template bottomRows<NX>(); }
		template<class Matrix> void set_xMax(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._ub.template bottomRows<NX>() = val; }
		decltype(auto) get_xMax(std::size_t i) const { return stage(i)._ub.template bottomRows<NX>(); }

		template<class Matrix> void set_uMin(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._lb.template topRows<NU>() = val; }
		decltype(auto) get_uMin(std::size_t i) const { return stage(i)._lb.template topRows<NU>(); }
		template<class Matrix> void set_uMax(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { stage(i)._ub.template topRows<NU>() = val; }
		decltype(auto) get_uMax(std::size_t i) const { return stage(i)._ub.template topRows<NU>(); }

		template<class Matrix> friend void set_zendMin(HPMPCProblem& p, Eigen::MatrixBase<Matrix> const& val) { p._zendMin = val; }
		friend StateVector const& get_zendMin(HPMPCProblem const& p) { return p._zendMin; }
		template<class Matrix> friend void set_zendMax(HPMPCProblem& p, Eigen::MatrixBase<Matrix> const& val) { p._zendMax = val; }
		friend StateVector const& get_zendMax(HPMPCProblem const& p) { return p._zendMax; }
		*/

	private:
		struct StageData
		{
			// Hessian = [R, S; S', Q]
			HPMPC_RMatrix _R;
			HPMPC_SMatrix _S;
			HPMPC_QMatrix _Q;

			// Gradient = [r; q]
			InputVector _r;
			StateVector _q;

			// Inter-stage equalities x_{k+1} = A x_k + B u_k + c_k
			HPMPC_AMatrix _A;
			HPMPC_BMatrix _B;
			StateVector _b;

			// Inequality constraints d_{min} <= C x_k + D u_k <= d_{max}
			HPMPC_CMatrix _C;
			HPMPC_DMatrix _D;
			StageConstraintVector _dMin;
			StageConstraintVector _dMax;

			// Bound constraints:
			// lb <= [u; x] <= ub
			StateInputVector _lb;
			StateInputVector _ub;
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

		// In the MPC case, the equality for stage 0 reads:
		// x[1] = A[0]*x[0] + B[0]*u[0] + c[0] = B[0]*u[0] + (A[0]*x[0] + c[0]) = B[0]*u[0] + b0,
		// b0 = A[0]*x[0].
		// and x[0] is not treated as an optimization variable by HPMPC.
		StateVector mutable _b0;

		// If x[0] is equality constrained and we eliminate it, then we need to
		// account for the x[0]' * S[0] * u[0] term:
		// x[0]' * S[0] * u[0] + r[0]' * u[0] =  (x[0]' * S[0] + r[0]') * u[0] = r0' * u[0],
		// where r0 = S[0]' * x[0] + r[0]
		InputVector mutable _r0;

		// 1 matrix of size NCT x NX.
		EndStageConstraintMatrix _C_end;

		// 1 vector of size NCT
		EndStageConstraintVector _d_end_min;

		// 1 vector of size NCT
		EndStageConstraintVector _d_end_max;

		// "A" data array for HPMPC
		std::vector<double const *> _A;

		// "B" data array for HPMPC
		std::vector<double const *> _B;

		// "b" data array for HPMPC
		std::vector<double const *> mutable _b;

		// "Q" data array for HPMPC
		std::vector<double const *> _Q;

		// "S" data array for HPMPC
		std::vector<double const *> _S;

		// "R" data array for HPMPC
		std::vector<double const *> _R;

		// "q" data array for HPMPC
		std::vector<double const *> _q;

		// "r" data array for HPMPC
		std::vector<double const *> mutable _r;

		// "lb" data array for HPMPC
		std::vector<double const *> _lb;

		// "ub" data array for HPMPC
		std::vector<double const *> _ub;

		// "C" data array for HPMPC
		std::vector<double const *> _C;

		// "D" data array for HPMPC
		std::vector<double const *> _D;

		// "lg" data array for HPMPC
		std::vector<double const *> _lg;

		// "ug" data array for HPMPC
		std::vector<double const *> _ug;
	};
}
