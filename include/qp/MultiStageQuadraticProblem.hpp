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
		static unsigned const NZ = NX + NU;
		static unsigned const NC = NC_;
		static unsigned const NCT = NCT_;

		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NU, 1> InputVector;
		typedef Eigen::Matrix<double, NZ, 1> StateInputVector;
		typedef Eigen::Matrix<double, NZ, NZ, Eigen::RowMajor> StageHessianMatrix;
		typedef Eigen::Matrix<double, NX, NX, Eigen::RowMajor> EndStageHessianMatrix;
		typedef Eigen::Matrix<double, NZ, 1> StageGradientVector;
		typedef Eigen::Matrix<double, NX, 1> EndStageGradientVector;
		typedef Eigen::Matrix<double, NX, NZ, Eigen::RowMajor> InterStageMatrix;
		typedef Eigen::Matrix<double, NX, 1> InterStageVector;
		typedef Eigen::Matrix<double, NC, NZ, Eigen::RowMajor> StageConstraintMatrix;
		typedef Eigen::Matrix<double, NC, 1> StageConstraintVector;
		typedef Eigen::Matrix<double, NCT, NX, Eigen::RowMajor> EndStageConstraintMatrix;
		typedef Eigen::Matrix<double, NCT, 1> EndStageConstraintVector;

		size_type nT() const { return _stage.size(); }
		static constexpr size_type nX() { return NX; }
		static constexpr size_type nZ() { return NZ; }
		static constexpr size_type nU() { return NU; }
		static constexpr size_type nD() { return NC; }
		static constexpr size_type nDT() { return NCT; }

		StageHessianMatrix& H(size_type i) { return stage(i)._H; }
		StageHessianMatrix const& H(size_type i) const { return stage(i)._H; }

		template <typename Matrix>
		friend void set_H(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& H) { p.H(i) = H; }

		EndStageHessianMatrix& Hend() { return _Hend; }
		EndStageHessianMatrix const& Hend() const {	return _Hend; }

		template <typename Matrix>
		friend void set_Hend(MultiStageQuadraticProblem& p, Eigen::MatrixBase<Matrix> const& Hend) { p.Hend() = Hend; }

		StageGradientVector& g(size_type i)	{ return stage(i)._g; }
		StageGradientVector const& g(size_type i) const	{ return stage(i)._g; }

		template <typename Matrix>
		friend void set_g(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& g) { p.g(i) = g; }

		EndStageGradientVector& gend() { return _gend; }
		EndStageGradientVector const& gend() const { return _gend;	}
		
		template <typename Matrix>
		friend void set_gend(MultiStageQuadraticProblem& p, Eigen::MatrixBase<Matrix> const& gend) { p.gend() = gend; }

		InterStageMatrix& C(size_type i) { return stage(i)._C; }
		InterStageMatrix const& C(size_type i) const { return stage(i)._C; }

		template <typename Matrix>
		friend void set_C(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& C) { p.C(i) = C; }

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
				
		InterStageVector& c(size_type i) { return stage(i)._c; }
		InterStageVector const& c(size_type i) const { return stage(i)._c; }

		template <typename Matrix>
		friend void set_c(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& c) { p.c(i) = c; }

		StateInputVector& zMin(size_type i)	{ return stage(i)._zMin; }
		StateInputVector const& zMin(size_type i) const { return stage(i)._zMin; }

		template <typename Matrix>
		friend void set_zMin(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& z_min) { p.zMin(i) = z_min; }

		StateVector& zendMin() { return _zendMin; }
		StateVector const& zendMin() const { return _zendMin; }
		
		template <typename Matrix>
		friend void set_zendMin(MultiStageQuadraticProblem& p, Eigen::MatrixBase<Matrix> const& val) { p.zendMin() = val; }

		StateInputVector& zMax(size_type i)	{ return stage(i)._zMax; }
		StateInputVector const& zMax(size_type i) const { return stage(i)._zMax; }
		
		template <typename Matrix>
		friend void set_zMax(MultiStageQuadraticProblem& p, std::size_t i, Eigen::MatrixBase<Matrix> const& z_max) { p.zMax(i) = z_max; }

		StateVector& zendMax() { return _zendMax; }
		StateVector const& zendMax() const { return _zendMax; }

		template <typename Matrix>
		friend void set_zendMax(MultiStageQuadraticProblem& p, Eigen::MatrixBase<Matrix> const& val) { p.zendMax() = val; }

		MultiStageQuadraticProblem(size_type nt) : _stage(nt) {}

	private:
		struct StageData
		{
			StageHessianMatrix _H;
			StageGradientVector _g;
			StageConstraintMatrix _D;
			StageConstraintVector _dMin;
			StageConstraintVector _dMax;
			InterStageMatrix _C;
			InterStageVector _c;
			StateInputVector _zMin;
			StateInputVector _zMax;
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
}

namespace camels
{
	using namespace tmpc;
}
