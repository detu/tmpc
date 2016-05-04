#pragma once

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
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class qpDUNESProgram
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

		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		void PrintQP_C(std::ostream& os) const;

		void PrintQP_zMax_C(std::ostream& log_stream) const;

		void PrintQP_zMin_C(std::ostream& log_stream) const;

		void PrintQP_MATLAB(std::ostream& log_stream, const std::string& var_name = "qp") const;

		size_type nT() const { return _Nt; }
		size_type nX() const { return NX; }
		size_type nZ() const { return NZ; }
		size_type nU() const { return NU; }
		size_type nD() const { return NC; }
		size_type nDT() const { return NCT; }

		size_type nIndep() const { return NX + NU * _Nt; }
		size_type nDep() const { return NX * _Nt; }
		size_type nVar() const { return NZ * _Nt + NX; }
		size_type nConstr() const { return NC * _Nt + NCT; }

		Eigen::Map<StageHessianMatrix> H(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StageHessianMatrix>(_H.data() + i * nZ() * nZ());
		}

		Eigen::Map<StageHessianMatrix const> H(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StageHessianMatrix const>(_H.data() + i * nZ() * nZ());
		}

		Eigen::Map<EndStageHessianMatrix> Hend()
		{
			return Eigen::Map<EndStageHessianMatrix>(_H.data() + nT() * nZ() * nZ());
		}

		Eigen::Map<EndStageHessianMatrix const> Hend() const
		{
			return Eigen::Map<EndStageHessianMatrix const>(_H.data() + nT() * nZ() * nZ());
		}

		Eigen::Map<StageGradientVector> g(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StageGradientVector>(_g.data() + i * nZ());
		}

		Eigen::Map<StageGradientVector const> g(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StageGradientVector const>(_g.data() + i * nZ());
		}

		Eigen::Map<EndStageGradientVector> gend()
		{
			return Eigen::Map<EndStageGradientVector>(_g.data() + nT() * nZ());
		}

		Eigen::Map<EndStageGradientVector const> gend() const
		{
			return Eigen::Map<EndStageGradientVector const>(_g.data() + nT() * nZ());
		}
		
		Eigen::Map<InterStageMatrix> C(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<InterStageMatrix>(_C.data() + i * nX() * nZ());
		}

		Eigen::Map<InterStageMatrix const> C(unsigned i) const
		{
			return Eigen::Map<InterStageMatrix const>(_C.data() + i * nX() * nZ());
		}

		Eigen::Map<StageConstraintMatrix> D(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintMatrix>(_D.data() + i * nD() * nZ());
		}

		Eigen::Map<StageConstraintMatrix const> D(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintMatrix const>(_D.data() + i * nD() * nZ());
		}

		Eigen::Map<EndStageConstraintMatrix> Dend()
		{
			return Eigen::Map<EndStageConstraintMatrix>(_D.data() + nT() * nD() * nZ());
		}

		Eigen::Map<EndStageConstraintMatrix const> Dend() const
		{
			return Eigen::Map<EndStageConstraintMatrix const>(_D.data() + nT() * nD() * nZ());
		}

		Eigen::Map<StageConstraintVector> dMin(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintVector>(_dMin.data() + i * nD());
		}

		Eigen::Map<StageConstraintVector const> dMin(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintVector const>(_dMin.data() + i * nD());
		}

		Eigen::Map<EndStageConstraintVector> dendMin()
		{
			return Eigen::Map<EndStageConstraintVector>(_dMin.data() + nT() * nD());
		}

		Eigen::Map<EndStageConstraintVector const> dendMin() const
		{
			return Eigen::Map<EndStageConstraintVector const>(_dMin.data() + nT() * nD());
		}

		Eigen::Map<StageConstraintVector> dMax(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintVector>(_dMax.data() + i * nD());
		}

		Eigen::Map<StageConstraintVector const> dMax(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StageConstraintVector const>(_dMax.data() + i * nD());
		}

		Eigen::Map<EndStageConstraintVector> dendMax()
		{
			return Eigen::Map<EndStageConstraintVector>(_dMax.data() + nT() * nD());
		}

		Eigen::Map<EndStageConstraintVector const> dendMax() const
		{
			return Eigen::Map<EndStageConstraintVector const>(_dMax.data() + nT() * nD());
		}
				
		Eigen::Map<InterStageVector> c(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<InterStageVector>(_c.data() + i * nX());
		}

		Eigen::Map<InterStageVector const> c(unsigned i) const
		{
			return Eigen::Map<InterStageVector const>(_c.data() + i * nX());
		}

		Eigen::Map<StateInputVector> zMin(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StateInputVector>(_zMin.data() + i * nZ());
		}

		Eigen::Map<StateInputVector const> zMin(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StateInputVector const>(_zMin.data() + i * nZ());
		}

		Eigen::Map<StateVector> zendMin()
		{
			return Eigen::Map<StateVector>(_zMin.data() + nT() * nZ());
		}

		Eigen::Map<StateVector const> zendMin() const
		{
			return Eigen::Map<StateVector const>(_zMin.data() + nT() * nZ());
		}
		
		Eigen::Map<StateInputVector> zMax(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<StateInputVector>(_zMax.data() + i * nZ());
		}

		Eigen::Map<StateInputVector const> zMax(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<StateInputVector const>(_zMax.data() + i * nZ());
		}
		
		Eigen::Map<StateVector> zendMax()
		{
			return Eigen::Map<StateVector>(_zMax.data() + nT() * nZ());
		}

		Eigen::Map<StateVector const> zendMax() const
		{
			return Eigen::Map<StateVector const>(_zMax.data() + nT() * nZ());
		}

		Eigen::Map<StateVector> xMin(unsigned i)
		{
			assert(i < nT() + 1);
			return Eigen::Map<StateVector>(_zMin.data() + i * nZ());
		}

		Eigen::Map<StateVector const> xMin(unsigned i) const
		{
			assert(i < nT() + 1);
			return Eigen::Map<StateVector const>(_zMin.data() + i * nZ());
		}

		Eigen::Map<StateVector> xMax(unsigned i)
		{
			assert(i < nT() + 1);
			return Eigen::Map<StateVector>(_zMax.data() + i * nZ());
		}

		Eigen::Map<StateVector const> xMax(unsigned i) const
		{
			assert(i < nT() + 1);
			return Eigen::Map<StateVector const>(_zMax.data() + i * nZ());
		}

		Eigen::Map<InputVector> uMin(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<InputVector>(_zMin.data() + i * nZ() + nX());
		}

		Eigen::Map<InputVector const> uMin(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<InputVector const>(_zMin.data() + i * nZ() + nX());
		}

		Eigen::Map<InputVector> uMax(unsigned i)
		{
			assert(i < nT());
			return Eigen::Map<InputVector>(_zMax.data() + i * nZ() + nX());
		}

		Eigen::Map<InputVector const> uMax(unsigned i) const
		{
			assert(i < nT());
			return Eigen::Map<InputVector const>(_zMax.data() + i * nZ() + nX());
		}

		qpDUNESProgram(size_type nt)
		:	_Nt(nt)
		,	_H(NZ * NZ * nt + NX * NX)
		,	_g(NZ * nt + NX)
		,	_C(NX * NZ * nt)
		,	_c(NX * nt)
		,	_D(NC * NZ * nt + NCT * NX)
		,	_dMin(NC * nt + NCT)
		,	_dMax(NC * nt + NCT)
		,	_zMin(NZ * nt + NX)
		,	_zMax(NZ * nt + NX)
		{
		}

		// Interface for accessing the data in qpDUNES format.
		//
		double const * H_data() const { return _H.data(); }
		double const * g_data() const { return _g.data(); }
		double const * C_data() const { return _C.data(); }
		double const * c_data() const { return _c.data(); }
		double const * D_data() const { return _D.data(); }
		double const * dMin_data() const { return _dMin.data(); }
		double const * dMax_data() const { return _dMax.data(); }
		double const * zMin_data() const { return _zMin.data(); }
		double const * zMax_data() const { return _zMax.data(); }

	private:
		// Private data members.
		//

		size_type const _Nt;

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
	};

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	inline void qpDUNESProgram<NX_, NU_, NC_, NCT_>::PrintQP_C(std::ostream& log_stream) const
	{
		using std::endl;

		Eigen::IOFormat C_format(Eigen::StreamPrecision, 0, ", ", ",\n", "", "", "", "");

		log_stream << "const double H[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << H(i).format(C_format) << "," << endl;
		log_stream << "};" << endl << endl;

		log_stream << "const double g[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << g(i).transpose().format(C_format) << "," << endl;
		log_stream << "};" << endl << endl;

		log_stream << "const double C[] = {" << endl;
		for (unsigned i = 0; i < nT(); ++i)
			log_stream << C(i).format(C_format) << "," << endl;
		log_stream << "};" << endl << endl;

		log_stream << "const double c[] = {" << endl;
		for (unsigned i = 0; i < nT(); ++i)
			log_stream << c(i).transpose().format(C_format) << "," << endl;
		log_stream << "};" << endl << endl;

		log_stream << "const double D[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << D(i).format(C_format) << "," << endl;
		log_stream << "};" << endl << endl;

		log_stream << "const double dMin[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << dMin(i).transpose().format(C_format) << ",";
		log_stream << "};" << endl << endl;

		log_stream << "const double dMax[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << dMax(i).transpose().format(C_format) << ",";
		log_stream << "};" << endl << endl;

		PrintQP_zMin_C(log_stream);
		PrintQP_zMax_C(log_stream);
	}

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	inline void qpDUNESProgram<NX_, NU_, NC_, NCT_>::PrintQP_MATLAB(std::ostream& log_stream, const std::string& var_name) const
	{
		using std::endl;

		for (unsigned k = 0; k <= nT(); ++k)
		{
			log_stream << var_name << ".H{" << k + 1 << "} = [..." << endl << H(k) << "];" << endl;
			log_stream << var_name << ".g{" << k + 1 << "} = [..." << endl << g(k) << "];" << endl;

			if (k < nT())
			{
				log_stream << var_name << ".C{" << k + 1 << "} = [..." << endl << C(k) << "];" << endl;
				log_stream << var_name << ".c{" << k + 1 << "} = [..." << endl << c(k) << "];" << endl;
			}

			log_stream << var_name << ".D{" << k + 1 << "} = [..." << endl << D(k) << "];" << endl;
			log_stream << var_name << ".dMin{" << k + 1 << "} = [..." << endl << dMin(k) << "];" << endl;
			log_stream << var_name << ".dMax{" << k + 1 << "} = [..." << endl << dMax(k) << "];" << endl;

			log_stream << var_name << ".zMin{" << k + 1 << "} = [..." << endl << zMin(k) << "];" << endl;
			log_stream << var_name << ".zMax{" << k + 1 << "} = [..." << endl << zMax(k) << "];" << endl;
		}
	}

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	inline void qpDUNESProgram<NX_, NU_, NC_, NCT_>::PrintQP_zMin_C(std::ostream& log_stream) const
	{
		using std::endl;

		Eigen::IOFormat C_format(Eigen::StreamPrecision, 0, ", ", ",\n", "", "", "", "");

		log_stream << "const double zLow[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << zMin(i).transpose().format(C_format) << "," << endl;
		log_stream << endl << "};" << endl << endl;
	}

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	inline void qpDUNESProgram<NX_, NU_, NC_, NCT_>::PrintQP_zMax_C(std::ostream& log_stream) const
	{
		using std::endl;

		Eigen::IOFormat C_format(Eigen::StreamPrecision, 0, ", ", ",\n", "", "", "", "");

		log_stream << "const double zUpp[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << zMax(i).transpose().format(C_format) << "," << endl;
		log_stream << endl << "};" << endl << endl;
	}
}
