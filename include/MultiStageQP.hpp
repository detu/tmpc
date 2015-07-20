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

		template<class TMatrix>
		void Calculate_M(Eigen::MatrixBase<TMatrix>& M) const
		{
			assert(M.rows() == nDep() && M.cols() == nIndep());
			M.setZero();

			for (unsigned i = 0; i < _Nt; ++i)
			{
				if (i == 0)
				{
					M.block(i * _Nx, 0, _Nx, _Nx + _Nu) = C(i);
				}
				else
				{
					M.block(i * _Nx, 0, _Nx, _Nx + i * _Nu + _Nu) << C(i).leftCols(_Nx) * M.block((i - 1) * _Nx, 0, _Nx, _Nx + i * _Nu), C(i).rightCols(_Nu);
				}
			}
		}

		template<class TMatrix>
		void Calculate_v(Eigen::MatrixBase<TMatrix>& v) const
		{
			assert(v.rows() == nDep() && v.cols() == 1);

			for (unsigned i = 0; i < _Nt; ++i)
			{
				if (i == 0)
				{
					v.middleRows(i * _Nx, _Nx) = c(i);
				}
				else
				{
					v.middleRows(i * _Nx, _Nx) = C(i).leftCols(_Nx) * v.middleRows((i - 1) * _Nx, _Nx) + c(i);
				}
			}
		}

		template<class Scalar>
		void Calculate_FullH(Eigen::SparseMatrix<Scalar>& H) const
		{
			assert(H.rows() == nVar() && H.cols() == nVar());

			typedef Eigen::Triplet<Scalar> Triplet;
			std::vector<Triplet> t;
			t.reserve(_Nt * _Nz * _Nz + _Nx * _Nx);

			for (unsigned k = 0; k <= _Nt; ++k)
			{
				const auto H_k = this->H(k);
				for (unsigned i = 0; i < H_k.rows(); ++i)
				{
					for (unsigned j = 0; j < H_k.cols(); ++j)
						t.emplace_back(i + k * _Nz, j + k * _Nz, H_k(i, j));
				}
			}

			H.setFromTriplets(t.begin(), t.end());
		}

		template<class Scalar>
		void Calculate_PermutationMatrix(Eigen::SparseMatrix<Scalar>& P) const
		{
			assert(P.rows() == _Nz * _Nt + _Nx && P.cols() == _Nz * _Nt + _Nx);
			P.setZero();
			P.reserve(_Nz * _Nt + _Nx);
			
			// x_0
			for (unsigned i = 0; i < _Nx; ++i)
				P.insert(i, i) = 1.;

			// x_k, k = 1...Nt
			for (unsigned k = 1; k <= _Nt; ++k)
				for (unsigned i = 0; i < _Nx; ++i)
					P.insert(k * _Nz + i, _Nt * _Nu + k * _Nx + i) = 1.;

			// u_k, k = 0...Nt-1
			for (unsigned k = 0; k < _Nt; ++k)
				for (unsigned i = 0; i < _Nu; ++i)
					P.insert(_Nx + k * _Nz + i, _Nx + k * _Nu + i) = 1.;
		}

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
		VectorMap xMax(unsigned i);

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
