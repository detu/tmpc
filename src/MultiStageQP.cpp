#include <MultiStageQP.hpp>

#include <stdexcept>

namespace camels
{
	MultiStageQP::MultiStageQP(size_type nx, size_type nu, size_type nd, size_type ndt, size_type nt) : 
		_Nx(nx),
		_Nu(nu),
		_Nt(nt),
		_Nz(nx + nu),
		_Nd(nd),
		_NdT(ndt)
	{
		// Allocate arrays.
		_H.resize(_Nz * _Nz * _Nt + _Nx * _Nx);
		_g.resize(_Nz * _Nt + _Nx);
		_C.resize(_Nx * _Nz * _Nt);
		_c.resize(_Nx * _Nt);
		_D.resize(_Nd * _Nz * _Nt + _NdT * _Nx);
		_dMin.resize(_Nd * _Nt + _NdT);
		_dMax.resize(_Nd * _Nt + _NdT);
		_zMin.resize(_Nz * _Nt + _Nx);
		_zMax.resize(_Nz * _Nt + _Nx);
		_zOpt.resize(_Nz * _Nt + _Nx);
	}

	MultiStageQP::~MultiStageQP()
	{
	}

	void MultiStageQP::PrintQP_C(std::ostream& log_stream) const
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

	void MultiStageQP::PrintQP_MATLAB(std::ostream& log_stream) const
	{
		using std::endl;

		for (unsigned k = 0; k <= _Nt; ++k)
		{
			log_stream << "qp.H{" << k + 1 << "} = [..." << endl << H(k) << "];" << endl;
			log_stream << "qp.g{" << k + 1 << "} = [..." << endl << g(k) << "];" << endl;
			
			if (k < _Nt)
			{
				log_stream << "qp.C{" << k + 1 << "} = [..." << endl << C(k) << "];" << endl;
				log_stream << "qp.c{" << k + 1 << "} = [..." << endl << c(k) << "];" << endl;
			}

			log_stream << "qp.D{" << k + 1 << "} = [..." << endl << D(k) << "];" << endl;
			log_stream << "qp.dMin{" << k + 1 << "} = [..." << endl << dMin(k) << "];" << endl;
			log_stream << "qp.dMax{" << k + 1 << "} = [..." << endl << dMax(k) << "];" << endl;

			log_stream << "qp.zMin{" << k + 1 << "} = [..." << endl << zMin(k) << "];" << endl;
			log_stream << "qp.zMax{" << k + 1 << "} = [..." << endl << zMax(k) << "];" << endl;
		}
	}

	MultiStageQP::RowMajorMatrixMap MultiStageQP::H(unsigned i)
	{
		assert(i < _Nt + 1);
		const auto sz = i < _Nt ? _Nz : _Nx;
		return RowMajorMatrixMap(_H.data() + i * _Nz * _Nz, sz, sz);
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::H(unsigned i) const
	{
		assert(i < _Nt + 1);
		const auto sz = i < _Nt ? _Nz : _Nx;
		return RowMajorMatrixConstMap(_H.data() + i * _Nz * _Nz, sz, sz);
	}

	MultiStageQP::VectorMap MultiStageQP::g(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_g.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::g(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_g.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::RowMajorMatrixMap MultiStageQP::C(unsigned i)
	{
		return RowMajorMatrixMap(_C.data() + i * _Nx * _Nz, _Nx, _Nz);
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::C(unsigned i) const
	{
		return RowMajorMatrixConstMap(_C.data() + i * _Nx * _Nz, _Nx, _Nz);
	}

	MultiStageQP::VectorMap MultiStageQP::c(unsigned i)
	{
		return VectorMap(_c.data() + i * _Nx, _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::c(unsigned i) const
	{
		return VectorConstMap(_c.data() + i * _Nx, _Nx);
	}

	MultiStageQP::VectorMap MultiStageQP::zMin(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMin.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::zMin(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMin.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::VectorMap MultiStageQP::zMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::zMax(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMax.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::VectorMap MultiStageQP::xMin(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMin.data() + i * _Nz, _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::xMin(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMin.data() + i * _Nz, _Nx);
	}

	MultiStageQP::VectorMap MultiStageQP::xMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, _Nx);
	}

	MultiStageQP::VectorConstMap MultiStageQP::xMax(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMax.data() + i * _Nz, _Nx);
	}

	void MultiStageQP::PrintQP_zMin_C(std::ostream& log_stream) const
	{
		using std::endl;

		Eigen::IOFormat C_format(Eigen::StreamPrecision, 0, ", ", ",\n", "", "", "", "");

		log_stream << "const double zLow[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << zMin(i).transpose().format(C_format) << "," << endl;
		log_stream << endl << "};" << endl << endl;
	}

	void MultiStageQP::PrintQP_zMax_C(std::ostream& log_stream) const
	{
		using std::endl;

		Eigen::IOFormat C_format(Eigen::StreamPrecision, 0, ", ", ",\n", "", "", "", "");

		log_stream << "const double zUpp[] = {" << endl;
		for (unsigned i = 0; i <= nT(); ++i)
			log_stream << zMax(i).transpose().format(C_format) << "," << endl;
		log_stream << endl << "};" << endl << endl;
	}

	MultiStageQP::VectorMap MultiStageQP::uMin(unsigned i)
	{
		assert(i < _Nt);
		return VectorMap(_zMin.data() + i * _Nz + _Nx, _Nu);
	}

	MultiStageQP::VectorConstMap MultiStageQP::uMin(unsigned i) const
	{
		assert(i < _Nt);
		return VectorConstMap(_zMin.data() + i * _Nz + _Nx, _Nu);
	}

	MultiStageQP::VectorMap MultiStageQP::uMax(unsigned i)
	{
		assert(i < _Nt);
		return VectorMap(_zMax.data() + i * _Nz + _Nx, _Nu);
	}

	MultiStageQP::VectorConstMap MultiStageQP::uMax(unsigned i) const
	{
		assert(i < _Nt);
		return VectorConstMap(_zMax.data() + i * _Nz + _Nx, _Nu);
	}

	MultiStageQP::VectorMap MultiStageQP::dMin(unsigned i)
	{
		assert(i <= _Nt);
		return VectorMap(_dMin.data() + i * _Nd, i < _Nt ? _Nd : _NdT);
	}

	MultiStageQP::VectorConstMap MultiStageQP::dMin(unsigned i) const
	{
		assert(i <= _Nt);
		return VectorConstMap(_dMin.data() + i * _Nd, i < _Nt ? _Nd : _NdT);
	}

	MultiStageQP::VectorMap MultiStageQP::dMax(unsigned i)
	{
		assert(i <= _Nt);
		return VectorMap(_dMax.data() + i * _Nd, i < _Nt ? _Nd : _NdT);
	}

	MultiStageQP::VectorConstMap MultiStageQP::dMax(unsigned i) const
	{
		assert(i <= _Nt);
		return VectorConstMap(_dMax.data() + i * _Nd, i < _Nt ? _Nd : _NdT);
	}

	MultiStageQP::RowMajorMatrixMap MultiStageQP::D(unsigned i)
	{
		assert(i <= _Nt);
		return RowMajorMatrixMap(_D.data() + i * _Nd * _Nz, i < _Nt ? _Nd : _NdT, i < _Nt ? _Nz : _Nx);
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::D(unsigned i) const
	{
		assert(i <= _Nt);
		return RowMajorMatrixConstMap(_D.data() + i * _Nd * _Nz, i < _Nt ? _Nd : _NdT, i < _Nt ? _Nz : _Nx);
	}
}
