#include <MultiStageQP.hpp>

#include <stdexcept>

namespace camels
{
	MultiStageQP::MultiStageQP(size_type nx, size_type nu, size_type nd, size_type ndt, size_type nt) :
		MultiStageQP(MultiStageQPSize(nx, nu, nd, ndt, nt))
	{
	}		

	MultiStageQP::MultiStageQP(const MultiStageQPSize& size) :
		_size(size)
	{
		// Allocate arrays.
		_H.resize(nZ() * nZ() * nT() + nX() * nX());
		_g.resize(nZ() * nT() + nX());
		_C.resize(nX() * nZ() * nT());
		_c.resize(nX() * nT());
		_D.resize(nD() * nZ() * nT() + nDT() * nX());
		_dMin.resize(nD() * nT() + nDT());
		_dMax.resize(nD() * nT() + nDT());
		_zMin.resize(nZ() * nT() + nX());
		_zMax.resize(nZ() * nT() + nX());
		_zOpt.resize(nZ() * nT() + nX());
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

	void MultiStageQP::PrintQP_MATLAB(std::ostream& log_stream, const std::string& var_name) const
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

	MultiStageQP::RowMajorMatrixMap MultiStageQP::H(unsigned i)
	{
		assert(i < nT() + 1);
		const auto sz = i < nT() ? nZ() : nX();
		return RowMajorMatrixMap(_H.data() + i * nZ() * nZ(), sz, sz);
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::H(unsigned i) const
	{
		assert(i < nT() + 1);
		const auto sz = i < nT() ? nZ() : nX();
		return RowMajorMatrixConstMap(_H.data() + i * nZ() * nZ(), sz, sz);
	}

	MultiStageQP::VectorMap MultiStageQP::g(unsigned i)
	{
		assert(i < nT() + 1);
		return VectorMap(_g.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::g(unsigned i) const
	{
		assert(i < nT() + 1);
		return VectorConstMap(_g.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::RowMajorMatrixMap MultiStageQP::C(unsigned i)
	{
		return RowMajorMatrixMap(_C.data() + i * nX() * nZ(), nX(), nZ());
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::C(unsigned i) const
	{
		return RowMajorMatrixConstMap(_C.data() + i * nX() * nZ(), nX(), nZ());
	}

	MultiStageQP::VectorMap MultiStageQP::c(unsigned i)
	{
		return VectorMap(_c.data() + i * nX(), nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::c(unsigned i) const
	{
		return VectorConstMap(_c.data() + i * nX(), nX());
	}

	MultiStageQP::VectorMap MultiStageQP::zMin(unsigned i)
	{
		assert(i < nT() + 1);
		return VectorMap(_zMin.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::zMin(unsigned i) const
	{
		assert(i < nT() + 1);
		return VectorConstMap(_zMin.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::VectorMap MultiStageQP::zMax(unsigned i)
	{
		assert(i < nT() + 1);
		return VectorMap(_zMax.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::zMax(unsigned i) const
	{
		assert(i < nT() + 1);
		return VectorConstMap(_zMax.data() + i * nZ(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::VectorMap MultiStageQP::xMin(unsigned i)
	{
		assert(i < nT() + 1);
		return VectorMap(_zMin.data() + i * nZ(), nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::xMin(unsigned i) const
	{
		assert(i < nT() + 1);
		return VectorConstMap(_zMin.data() + i * nZ(), nX());
	}

	MultiStageQP::VectorMap MultiStageQP::xMax(unsigned i)
	{
		assert(i < nT() + 1);
		return VectorMap(_zMax.data() + i * nZ(), nX());
	}

	MultiStageQP::VectorConstMap MultiStageQP::xMax(unsigned i) const
	{
		assert(i < nT() + 1);
		return VectorConstMap(_zMax.data() + i * nZ(), nX());
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
		assert(i < nT());
		return VectorMap(_zMin.data() + i * nZ() + nX(), nU());
	}

	MultiStageQP::VectorConstMap MultiStageQP::uMin(unsigned i) const
	{
		assert(i < nT());
		return VectorConstMap(_zMin.data() + i * nZ() + nX(), nU());
	}

	MultiStageQP::VectorMap MultiStageQP::uMax(unsigned i)
	{
		assert(i < nT());
		return VectorMap(_zMax.data() + i * nZ() + nX(), nU());
	}

	MultiStageQP::VectorConstMap MultiStageQP::uMax(unsigned i) const
	{
		assert(i < nT());
		return VectorConstMap(_zMax.data() + i * nZ() + nX(), nU());
	}

	MultiStageQP::VectorMap MultiStageQP::dMin(unsigned i)
	{
		assert(i <= nT());
		return VectorMap(_dMin.data() + i * nD(), i < nT() ? nD() : nDT());
	}

	MultiStageQP::VectorConstMap MultiStageQP::dMin(unsigned i) const
	{
		assert(i <= nT());
		return VectorConstMap(_dMin.data() + i * nD(), i < nT() ? nD() : nDT());
	}

	MultiStageQP::VectorMap MultiStageQP::dMax(unsigned i)
	{
		assert(i <= nT());
		return VectorMap(_dMax.data() + i * nD(), i < nT() ? nD() : nDT());
	}

	MultiStageQP::VectorConstMap MultiStageQP::dMax(unsigned i) const
	{
		assert(i <= nT());
		return VectorConstMap(_dMax.data() + i * nD(), i < nT() ? nD() : nDT());
	}

	MultiStageQP::RowMajorMatrixMap MultiStageQP::D(unsigned i)
	{
		assert(i <= nT());
		return RowMajorMatrixMap(_D.data() + i * nD() * nZ(), i < nT() ? nD() : nDT(), i < nT() ? nZ() : nX());
	}

	MultiStageQP::RowMajorMatrixConstMap MultiStageQP::D(unsigned i) const
	{
		assert(i <= nT());
		return RowMajorMatrixConstMap(_D.data() + i * nD() * nZ(), i < nT() ? nD() : nDT(), i < nT() ? nZ() : nX());
	}

	camels::MultiStageQP::size_type MultiStageQP::nT() const
	{
		return _size.nT();
	}

	camels::MultiStageQP::size_type MultiStageQP::nX() const
	{
		return _size.nX();
	}

	camels::MultiStageQP::size_type MultiStageQP::nZ() const
	{
		return _size.nZ();
	}

	camels::MultiStageQP::size_type MultiStageQP::nU() const
	{
		return _size.nU();
	}

	camels::MultiStageQP::size_type MultiStageQP::nD() const
	{
		return _size.nD();
	}

	camels::MultiStageQP::size_type MultiStageQP::nDT() const
	{
		return _size.nDT();
	}

	camels::MultiStageQP::size_type MultiStageQP::nIndep() const
	{
		return _size.nIndep();
	}

	camels::MultiStageQP::size_type MultiStageQP::nDep() const
	{
		return _size.nDep();
	}

	camels::MultiStageQP::size_type MultiStageQP::nVar() const
	{
		return _size.nVar();
	}

	const MultiStageQPSize& MultiStageQP::size() const
	{
		return _size;
	}
}
