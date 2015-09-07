#pragma once

#include "MultiStageQPSize.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <vector>
#include <memory>
#include <ostream>
#include <cstdlib>

namespace camels
{
	/* MultiStageQP represents a problem
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
	class MultiStageQP
	{
	public:
		typedef unsigned int size_type;

		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		MultiStageQP(const MultiStageQPSize& size);
		MultiStageQP(size_type nx, size_type nu, size_type nd, size_type ndt, size_type nt);
		~MultiStageQP();

		void PrintQP_C(std::ostream& os) const;

		void PrintQP_zMax_C(std::ostream& log_stream) const;

		void PrintQP_zMin_C(std::ostream& log_stream) const;

		void PrintQP_MATLAB(std::ostream& log_stream) const;

		const MultiStageQPSize& size() const;
		size_type nT() const;
		size_type nX() const;
		size_type nZ() const;
		size_type nU() const;
		size_type nD() const;
		size_type nDT() const;
		size_type nIndep() const;
		size_type nDep() const;
		size_type nVar() const;

		RowMajorMatrixMap H(unsigned i);
		RowMajorMatrixConstMap H(unsigned i) const;
		
		VectorMap g(unsigned i);
		VectorConstMap g(unsigned i) const;
		
		RowMajorMatrixMap C(unsigned i);
		RowMajorMatrixConstMap C(unsigned i) const;

		RowMajorMatrixMap D(unsigned i);
		RowMajorMatrixConstMap D(unsigned i) const;

		VectorMap dMin(unsigned i);
		VectorConstMap dMin(unsigned i) const;

		VectorMap dMax(unsigned i);
		VectorConstMap dMax(unsigned i) const;
				
		VectorMap c(unsigned i);
		VectorConstMap c(unsigned i) const;
		
		VectorMap zMin(unsigned i);
		VectorConstMap zMin(unsigned i) const;
		
		VectorMap zMax(unsigned i);
		VectorConstMap zMax(unsigned i) const;

		VectorMap xMin(unsigned i);
		VectorConstMap xMin(unsigned i) const;

		VectorMap xMax(unsigned i);
		VectorConstMap xMax(unsigned i) const;

		VectorMap uMin(unsigned i);
		VectorConstMap uMin(unsigned i) const;

		VectorMap uMax(unsigned i);
		VectorConstMap uMax(unsigned i) const;

	private:
		const MultiStageQPSize _size;

		// _H stores _Nt row-major matrices of size _Nz x _Nz and 1 matrix of size _Nx x _Nx.
		std::vector<double> _H;

		// _g stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _g;

		// _C stores _Nt row-major matrices of size _Nx x _Nz
		std::vector<double> _C;

		// _c stores _Nt vectors of size _Nx
		std::vector<double> _c;

		// _D stores _Nt row-major matrices of size _Nd x _Nz and 1 row-major matrix of size _NdT x _Nx
		std::vector<double> _D;

		// _dMin stores _Nt vectors of size _Nd and 1 vector of size _NdT
		std::vector<double> _dMin;

		// _dMax stores _Nt vectors of size _Nd and 1 vector of size _NdT
		std::vector<double> _dMax;

		// _zMin stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMin;

		// _zMax stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMax;

		// Primal optimal solution.
		// _zOpt stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zOpt;
	};
}
