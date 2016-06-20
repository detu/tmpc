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

		size_type nT() const { return _stage.size(); }
		static constexpr size_type nX() { return NX; }
		static constexpr size_type nZ() { return NZ; }
		static constexpr size_type nU() { return NU; }
		static constexpr size_type nD() { return NC; }
		static constexpr size_type nDT() { return NCT; }

		template<class Matrix> void set_Q(std::size_t i, Eigen::MatrixBase<Matrix> const& Q) { stage(i)._Q = Q; }
		HPMPC_QMatrix const& get_Q(std::size_t i) const { return stage(i)._Q; }

		template<class Matrix> void set_R(std::size_t i, Eigen::MatrixBase<Matrix> const& R) { stage(i)._R = R; }
		HPMPC_RMatrix const& get_R(std::size_t i) const { return stage(i)._R; }

		template<class Matrix> void set_S(std::size_t i, Eigen::MatrixBase<Matrix> const& S) { stage(i)._S = S; }
		HPMPC_SMatrix const& get_S(std::size_t i) const { return stage(i)._S; }

		template<class Matrix> void set_A(std::size_t i, Eigen::MatrixBase<Matrix> const& A) { stage(i)._A = A; }
		HPMPC_AMatrix const& get_A(std::size_t i) const { return stage(i)._A; }

		template<class Matrix> void set_B(std::size_t i, Eigen::MatrixBase<Matrix> const& B) { stage(i)._B = B; }
		HPMPC_BMatrix const& get_B(std::size_t i) const { return stage(i)._B; }

		template<class Matrix> friend void set_Hend(HPMPCProblem& p, Eigen::MatrixBase<Matrix> const& Hend) { p._Hend = Hend; }
		friend EndStageHessianMatrix const& get_Hend(HPMPCProblem const& p) { return p._Hend; }

		template<class Matrix> friend void set_g(HPMPCProblem& p, size_type i, Eigen::MatrixBase<Matrix>& val)	{ p.stage(i)._g = val; }
		friend StageGradientVector const& get_g(HPMPCProblem const& p, size_type i) { return p.stage(i)._g; }

		template<class Matrix> friend void set_gend(HPMPCProblem& p, Eigen::MatrixBase<Matrix>& val) { p._gend = val; }
		friend EndStageGradientVector const& get_gend(HPMPCProblem const& p) { return p._gend;	}

		/*
		StageConstraintMatrix& D(size_type i) {	return stage(i)._D; }
		StageConstraintMatrix const& D(size_type i) const {	return stage(i)._D; }

		EndStageConstraintMatrix& Dend() { return _Dend; }
		EndStageConstraintMatrix const& Dend() const { return _Dend; }

		StageConstraintVector& dMin(size_type i) { return stage(i)._dMin; }
		StageConstraintVector const& dMin(size_type i) const { return stage(i)._dMin; }

		EndStageConstraintVector& dendMin() { return _dendMin; }
		EndStageConstraintVector const& dendMin() const	{ return _dendMin; }

		StageConstraintVector& dMax(size_type i) { return stage(i)._dMax; }
		StageConstraintVector const& dMax(size_type i) const { return stage(i)._dMax; }

		EndStageConstraintVector& dendMax()	{ return _dendMax; }
		EndStageConstraintVector const& dendMax() const	{ return _dendMax; }
		*/

		template<class Matrix> friend void set_c(HPMPCProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& c) { p.stage(i)._c = c; }
		friend StateVector const& get_c(HPMPCProblem& p, std::size_t i) { return p.stage(i)._c; }

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

		// ******************************************************
		//                HPMPC raw data interface.
		//
		// The prefixes before _data() correspond to the names of
		// the argument to c_order_d_ip_ocp_hard_tv().
		// ******************************************************
		double const * const * A_data () const { return _A .data(); }
		double const * const * B_data () const { return _B .data(); }
		double const * const * b_data () const { return _b .data(); }
		double const * const * Q_data () const { return _Q .data(); }
		double const * const * S_data () const { return _S .data(); }
		double const * const * R_data () const { return _R .data(); }
		double const * const * q_data () const { return _q .data(); }
		double const * const * r_data () const { return _r .data(); }
		double const * const * lb_data() const { return _lb.data(); }
		double const * const * ub_data() const { return _ub.data(); }
		double const * const * C_data () const { return _C .data(); }
		double const * const * D_data () const { return _D .data(); }
		double const * const * lg_data() const { return _lg.data(); }
		double const * const * ug_data() const { return _ug.data(); }

		HPMPCProblem(size_type nt)
		:	_stage(nt)
		,	_A(nt)
		,	_B(nt)
		,	_b(nt)
		,	_Q(nt + 1)
		,	_S(nt + 1)
		,	_R(nt + 1)
		,	_q(nt + 1)
		,	_r(nt + 1)
		,	_lb(nt + 1)
		,	_ub(nt + 1)
		,	_C(nt + 1)
		,	_D(nt + 1)
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
					_b[i] = _stage[i]._c.data();
				}

				_Q [i] = i < nt ? _stage[i]._Q.data() : _Hend.data();
				_S [i] = i < nt ? _stage[i]._S .data() : nullptr;
				_R [i] = i < nt ? _stage[i]._R .data() : nullptr;
				_q [i] = i < nt ? _stage[i]._q .data() : _gend.data();
				_r [i] = i < nt ? _stage[i]._r .data() : _gend.data() + NX;
				_lb[i] = i < nt ? _stage[i]._lb.data() : _zendMin.data();
				_ub[i] = i < nt ? _stage[i]._ub.data() : _zendMax.data();
				_C [i] = i < nt ? _stage[i]._C .data() : _Dend.data();
				_D [i] = i < nt ? _stage[i]._D .data() : nullptr;
				_lg[i] = i < nt ? _stage[i]._dMin.data() : _dendMin.data();
				_ug[i] = i < nt ? _stage[i]._dMax.data() : _dendMax.data();
			}
		}

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
			StateVector _c;

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

		StageData& stage(size_type i)
		{
			assert(i < nT());
			return _stage[i];
		}

		StageData const& stage(size_type i) const
		{
			assert(i < nT());
			return _stage[i];
		}

		// Private data members.
		//

		// Stores stage data
		std::vector<StageData> _stage;

		// 1 matrix of size _Nx x _Nx.
		EndStageHessianMatrix _Hend;

		// 1 vector of size _Nx
		EndStageGradientVector _gend;

		// 1 matrix of size NCT x NX.
		EndStageConstraintMatrix _Dend;

		// 1 vector of size NCT
		EndStageConstraintVector _dendMin;

		// 1 vector of size NCT
		EndStageConstraintVector _dendMax;

		// 1 vector of size _Nx
		StateVector _zendMin;

		// 1 vector of size _Nx
		StateVector _zendMax;

		// "A" data array for HPMPC
		std::vector<double const *> _A;

		// "B" data array for HPMPC
		std::vector<double const *> _B;

		// "b" data array for HPMPC
		std::vector<double const *> _b;

		// "Q" data array for HPMPC
		std::vector<double const *> _Q;

		// "S" data array for HPMPC
		std::vector<double const *> _S;

		// "R" data array for HPMPC
		std::vector<double const *> _R;

		// "q" data array for HPMPC
		std::vector<double const *> _q;

		// "r" data array for HPMPC
		std::vector<double const *> _r;

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

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void set_H(HPMPCProblem<NX, NU, NC, NCT>& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == NX + NU && Matrix::ColsAtCompileTime == NX + NU,
			"Matrix of size (NX+NU)x(NX+NU) is expected");
		qp.set_Q(i, val.template topLeftCorner<NX, NX>());
		qp.set_S(i, val.template bottomLeftCorner<NU, NX>());
		qp.set_R(i, val.template bottomRightCorner<NU, NU>());
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT>
	typename HPMPCProblem<NX, NU, NC, NCT>::StageHessianMatrix get_H(HPMPCProblem<NX, NU, NC, NCT> const& qp, std::size_t i)
	{
		typename HPMPCProblem<NX, NU, NC, NCT>::StageHessianMatrix H;
		H << qp.get_Q(i), qp.get_S(i).transpose(), qp.get_S(i), qp.get_R(i);
		return H;
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void set_C(HPMPCProblem<NX, NU, NC, NCT>& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& C)
	{
		static_assert(Matrix::ColsAtCompileTime == NX + NU,	"Matrix with (NX+NU) columns is expected");
		qp.set_A(i, C.template leftCols<NX>());
		qp.set_B(i, C.template rightCols<NU>());
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT>
	typename HPMPCProblem<NX, NU, NC, NCT>::InterStageMatrix get_C(HPMPCProblem<NX, NU, NC, NCT> const& qp, std::size_t i)
	{
		typename HPMPCProblem<NX, NU, NC, NCT>::InterStageMatrix C;
		C << qp.get_A(i), qp.get_B(i);
		return C;
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void set_zMin(HPMPCProblem<NX, NU, NC, NCT>& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& z_min)
	{
		static_assert(Matrix::RowsAtCompileTime == NX + NU,	"Vector with (NX+NU) rows is expected");
		qp.set_xMin(i, z_min.template topRows<NX>());
		qp.set_uMin(i, z_min.template bottomRows<NU>());
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT>
	typename HPMPCProblem<NX, NU, NC, NCT>::StateInputVector get_zMin(HPMPCProblem<NX, NU, NC, NCT> const& qp, std::size_t i)
	{
		typename HPMPCProblem<NX, NU, NC, NCT>::StateInputVector z_min;
		z_min << qp.get_xMin(i), qp.get_uMin(i);
		return z_min;
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void set_zMax(HPMPCProblem<NX, NU, NC, NCT>& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& z_max)
	{
		static_assert(Matrix::RowsAtCompileTime == NX + NU,	"Vector with (NX+NU) rows is expected");
		qp.set_xMax(i, z_max.template topRows<NX>());
		qp.set_uMax(i, z_max.template bottomRows<NU>());
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT>
	typename HPMPCProblem<NX, NU, NC, NCT>::StateInputVector get_zMax(HPMPCProblem<NX, NU, NC, NCT> const& qp, std::size_t i)
	{
		typename HPMPCProblem<NX, NU, NC, NCT>::StateInputVector z_max;
		z_max << qp.get_xMax(i), qp.get_uMax(i);
		return z_max;
	}
}
