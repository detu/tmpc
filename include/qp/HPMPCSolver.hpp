#pragma once

#include "MultiStageQuadraticProblemBase.hpp"
#include "qpDUNESSolution.hpp"

#include <c_interface.h>	// TODO: "c_interface.h" is a very general name; change the design so that it is <hpmpc/c_interface.h>, for example.

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

		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NU, 1> InputVector;
		typedef Eigen::Matrix<double, NZ, 1> StateInputVector;
		typedef Eigen::Matrix<double, NZ, NZ, Eigen::ColMajor> StageHessianMatrix;
		typedef Eigen::Matrix<double, NX, NX, Eigen::ColMajor> EndStageHessianMatrix;
		typedef Eigen::Matrix<double, NZ, 1> StageGradientVector;
		typedef Eigen::Matrix<double, NX, 1> EndStageGradientVector;
		typedef Eigen::Matrix<double, NX, NZ, Eigen::ColMajor> InterStageMatrix;
		typedef Eigen::Matrix<double, NC, NZ, Eigen::ColMajor> StageConstraintMatrix;
		typedef Eigen::Matrix<double, NC, 1> StageConstraintVector;
		typedef Eigen::Matrix<double, NCT, NX, Eigen::ColMajor> EndStageConstraintMatrix;
		typedef Eigen::Matrix<double, NCT, 1> EndStageConstraintVector;

		size_type nT() const { return _stage.size(); }
		static constexpr size_type nX() { return NX; }
		static constexpr size_type nZ() { return NZ; }
		static constexpr size_type nU() { return NU; }
		static constexpr size_type nD() { return NC; }
		static constexpr size_type nDT() { return NCT; }

		StageHessianMatrix& H(size_type i) { return stage(i)._H; }
		StageHessianMatrix const& H(size_type i) const { return stage(i)._H; }

		EndStageHessianMatrix& Hend() { return _Hend; }
		EndStageHessianMatrix const& Hend() const {	return _Hend; }

		StageGradientVector& g(size_type i)	{ return stage(i)._g; }
		StageGradientVector const& g(size_type i) const	{ return stage(i)._g; }

		EndStageGradientVector& gend() { return _gend; }
		EndStageGradientVector const& gend() const { return _gend;	}

		InterStageMatrix& C(size_type i) { return stage(i)._C; }
		InterStageMatrix const& C(size_type i) const { return stage(i)._C; }

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

		StateVector& c(size_type i) { return stage(i)._c; }
		StateVector const& c(size_type i) const { return stage(i)._c; }

		friend void setXMin(HPMPCProblem& p, std::size_t i, StateVector const& val) { p.stage(i)._xMin = val; }
		friend void setUMin(HPMPCProblem& p, std::size_t i, InputVector const& val) { p.stage(i)._uMin = val; }
		StateInputVector const& zMin(size_type i) const { return stage(i)._zMin; }

		friend void setZEndMin(HPMPCProblem& p, StateVector const& val) { p._zendMin = val; }
		StateVector const& zendMin() const { return _zendMin; }

		friend void setXMax(HPMPCProblem& p, std::size_t i, StateVector const& val) { p.stage(i)._xMax = val; }
		friend void setUMax(HPMPCProblem& p, std::size_t i, InputVector const& val) { p.stage(i)._uMax = val; }
		StateInputVector const& zMax(size_type i) const { return stage(i)._zMax; }

		friend void setZEndMax(HPMPCProblem& p, StateVector const& val) { p._zendMax = val; }
		StateVector const& zendMax() const { return _zendMax; }

		HPMPCProblem(size_type nt) : _stage(nt) {}

	private:
		typedef HPMPCMatrix<NU, NU> HPMPC_RMatrix;
		typedef HPMPCMatrix<NU, NX> HPMPC_SMatrix;
		typedef HPMPCMatrix<NX, NX> HPMPC_QMatrix;
		typedef HPMPCMatrix<NX, NX> HPMPC_AMatrix;
		typedef HPMPCMatrix<NX, NU> HPMPC_BMatrix;
		typedef HPMPCMatrix<NC, NX> HPMPC_CMatrix;
		typedef HPMPCMatrix<NC, NU> HPMPC_DMatrix;

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

			// Bound constraints
			// u_{min} <= u <= u_{max}
			// x_{min} <= x <= x_{max}
			InputVector _uMin;
			InputVector _uMax;
			StateVector _xMin;
			StateVector _xMax;
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
	};

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class HPMPCSolver
	{
	public:
		typedef HPMPCProblem<NX_, NU_, NC_, NCT_> Problem;
		typedef qpDUNESSolution<NX_, NU_> Solution;
	};

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void setZMin(HPMPCProblem<NX, NU, NC, NCT>& qp, size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		setXMin(qp, i, val.template topRows<NX>());
		setUMin(qp, i, val.template bottomRows<NU>());
	}

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT, typename Matrix>
	void setZMax(HPMPCProblem<NX, NU, NC, NCT>& qp, size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		setXMax(qp, i, val.template topRows<NX>());
		setUMax(qp, i, val.template bottomRows<NU>());
	}
}
