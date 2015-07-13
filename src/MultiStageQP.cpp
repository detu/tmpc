#include <MultiStageQP.hpp>

#include <stdexcept>

namespace camels
{
	MultiStageQP::MultiStageQP(size_type nx, size_type nu, size_type nt)
		: _levenbergMarquardt(0.01)
		, _Nx(nx)
		, _Nu(nu)
		, _Nt(nt)
		, _Nz(nx + nu)
	{
		// Allocate arrays.
		_H.resize(_Nz * _Nz * _Nt + _Nx * _Nx);
		_g.resize(_Nz * _Nt + _Nx);
		_C.resize(_Nx * _Nz * _Nt);
		_c.resize(_Nx * _Nt);
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

		log_stream << "const double H[] = {" << endl;
		for (unsigned i = 0; i < _H.size(); ++i)
		{
			log_stream << _H[i] << ",\t";

			if (i < _Nz * _Nz * _Nt)
			{
				if ((i + 1) % _Nz == 0)
					log_stream << endl;

				if ((i + 1) % (_Nz * _Nz) == 0)
					log_stream << endl;
			}
			else
			{
				if ((i - _Nz * _Nz * _Nt + 1) % _Nx == 0)
					log_stream << endl;

				if ((i - _Nz * _Nz * _Nt + 1) % (_Nx * _Nx) == 0)
					log_stream << endl;
			}
		}
		log_stream << "};" << endl << endl;

		log_stream << "const double g[] = {" << endl;
		for (unsigned i = 0; i < _g.size(); ++i)
		{
			log_stream << _g[i] << ",\t";

			if (i < _Nz * _Nt)
			{
				if ((i + 1) % _Nz == 0)
					log_stream << endl;
			}
			else
			{
				if ((i - _Nz * _Nt + 1) % _Nx == 0)
					log_stream << endl;
			}
		}
		log_stream << "};" << endl << endl;

		log_stream << "const double C[] = {" << endl;
		for (unsigned i = 0; i < _C.size(); ++i)
		{
			log_stream << _C[i] << ",\t";

			if ((i + 1) % _Nz == 0)
				log_stream << endl;

			if ((i + 1) % (_Nz * _Nx) == 0)
				log_stream << endl;
		}
		log_stream << "};" << endl << endl;

		log_stream << "const double c[] = {" << endl;
		for (unsigned i = 0; i < _c.size(); ++i)
		{
			log_stream << _c[i] << ",\t";

			if ((i + 1) % _Nx == 0)
				log_stream << endl;
		}
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

	MultiStageQP::VectorMap MultiStageQP::xMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, _Nx);
	}

	void MultiStageQP::PrintQP_zMin_C(std::ostream& log_stream) const
	{
		using std::endl;

		log_stream << "const double zLow[] = {" << endl;
		for (unsigned i = 0; i < _zMin.size(); ++i)
		{
			log_stream << _zMin[i] << ",\t";

			if ((i + 1) % _Nz == 0)
				log_stream << endl;
		}
		log_stream << endl << "};" << endl << endl;
	}

	void MultiStageQP::PrintQP_zMax_C(std::ostream& log_stream) const
	{
		using std::endl;

		log_stream << "const double zUpp[] = {" << endl;
		for (unsigned i = 0; i < _zMax.size(); ++i)
		{
			log_stream << _zMax[i] << ",\t";

			if ((i + 1) % _Nz == 0)
				log_stream << endl;
		}
		log_stream << endl << "};" << endl << endl;
	}
}
