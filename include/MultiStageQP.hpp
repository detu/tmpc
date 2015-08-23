#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <vector>
#include <memory>
#include <ostream>
#include <cstdlib>

namespace camels
{
	class MultiStageQP
	{
	public:
		typedef unsigned int size_type;

		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		MultiStageQP(size_type nx, size_type nu, size_type nt);
		~MultiStageQP();

		void PrintQP_C(std::ostream& os) const;

		void PrintQP_zMax_C(std::ostream& log_stream) const;

		void PrintQP_zMin_C(std::ostream& log_stream) const;

		void PrintQP_MATLAB(std::ostream& log_stream) const;
		size_type nT() const { return _Nt; }
		size_type nX() const { return _Nx; }
		size_type nZ() const { return _Nz; }
		size_type nU() const { return _Nu; }
		size_type nIndep() const { return _Nx + _Nu * _Nt; }
		size_type nDep() const { return _Nx * _Nt; }
		size_type nVar() const { return _Nz * _Nt + _Nx; }

		RowMajorMatrixMap H(unsigned i);
		RowMajorMatrixConstMap H(unsigned i) const;
		
		VectorMap g(unsigned i);
		VectorConstMap g(unsigned i) const;
		
		RowMajorMatrixMap C(unsigned i);
		RowMajorMatrixConstMap C(unsigned i) const;

// 		RowMajorMatrixMap A(size_type i);
// 		RowMajorMatrixConstMap A(size_type i) const;
		
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
		size_type _Nu;
		size_type _Nx;
		size_type _Nz;
		size_type _Nt;

		double _levenbergMarquardt;
		
		// _H stores _Nt row-major matrices of size _Nz x _Nz and 1 matrix of size _Nx x _Nx.
		std::vector<double> _H;

		// _g stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _g;

		// _C stores _Nt row-major matrices of size _Nx x _Nz
		std::vector<double> _C;

		// _c stores _Nt vectors of size _Nx
		std::vector<double> _c;

		// _zMin stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMin;

		// _zMax stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMax;

		// Primal optimal solution.
		// _zOpt stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zOpt;
	};
}
