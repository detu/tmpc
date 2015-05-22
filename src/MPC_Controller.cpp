#include <MPC_Controller.hpp>

#include <stdexcept>

namespace rtmc
{
	Eigen::MatrixXd OutputWeightingMatrix()
	{
		Eigen::MatrixXd w(6, 6);
		w.fill(0.);
		w.diagonal()[0] = 1;
		w.diagonal()[1] = 1;
		w.diagonal()[2] = 1;
		w.diagonal()[3] = 10;
		w.diagonal()[4] = 10;
		w.diagonal()[5] = 10;

		return w.transpose() * w;
	}

	MPC_Controller::MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, 
		unsigned Nt, const qpOptions_t * qpOptions /*= nullptr*/) : _platform(platform)
		, _levenbergMarquardt(0.01)
		, _sampleTime(sample_time)
	{
		// Get sizes.
		_Nq = platform->getNumberOfAxes();
		_Nu = platform->getInputDim();
		_Nx = platform->getStateDim();
		_Ny = platform->getOutputDim();
		_Nz = _Nx + _Nu;
		_Nt = Nt;

		// Allocate arrays.		
		_H.resize(_Nz * _Nz * _Nt + _Nx * _Nx);
		_g.resize(_Nz * _Nt + _Nx);
		_C.resize(_Nx * _Nz * _Nt);
		_c.resize(_Nx * _Nt);
		_zMin.resize(_Nz * _Nt + _Nx);
		_zMax.resize(_Nz * _Nt + _Nx);
		_z.resize(_Nz * _Nt + _Nx);
		_yRef.resize(_Ny * _Nt);

		// Set up qpData
		unsigned int* nD = 0;	  			/* number of affine constraints */
		return_t statusFlag = qpDUNES_setup(&_qpData, _Nt, _Nx, _Nu, nD, qpOptions);
		
		if (statusFlag != QPDUNES_OK)
			throw std::runtime_error("qpDUNES_setup() failed.");

		InitWorkingPoint();

	}

	MPC_Controller::~MPC_Controller()
	{
		/** cleanup of allocated data */
		qpDUNES_cleanup(&_qpData);
	}

	void MPC_Controller::PrintQP(std::ostream& log_stream) const
	{
		log_stream << "H = " << std::endl;
		for (unsigned i = 0; i < _H.size(); ++i)
		{
			log_stream << _H[i] << '\t';

			if (i < _Nz * _Nz * _Nt)
			{
				if ((i + 1) % _Nz == 0)
					log_stream << std::endl;

				if ((i + 1) % (_Nz * _Nz) == 0)
					log_stream << std::endl;
			}
			else
			{
				if ((i - _Nz * _Nz * _Nt + 1) % _Nx == 0)
					log_stream << std::endl;

				if ((i - _Nz * _Nz * _Nt + 1) % (_Nx * _Nx) == 0)
					log_stream << std::endl;
			}
		}

		log_stream << "g = " << std::endl;
		for (unsigned i = 0; i < _g.size(); ++i)
		{
			log_stream << _g[i] << '\t';

			if (i < _Nz * _Nt)
			{
				if ((i + 1) % _Nz == 0)
					log_stream << std::endl;
			}
			else
			{
				if ((i - _Nz * _Nt + 1) % _Nx == 0)
					log_stream << std::endl;
			}
		}
		log_stream << std::endl;

		log_stream << "C = " << std::endl;
		for (unsigned i = 0; i < _C.size(); ++i)
		{
			log_stream << _C[i] << '\t';

			if ((i + 1) % _Nz == 0)
				log_stream << std::endl;

			if ((i + 1) % (_Nz * _Nx) == 0)
				log_stream << std::endl;
		}

		log_stream << "c = " << std::endl;
		for (unsigned i = 0; i < _c.size(); ++i)
		{
			log_stream << _c[i] << '\t';

			if ((i + 1) % _Nx == 0)
				log_stream << std::endl;
		}
		log_stream << std::endl;

		log_stream << "_zMin = " << std::endl;
		for (unsigned i = 0; i < _zMin.size(); ++i)
		{
			log_stream << _zMin[i] << '\t';

			if ((i + 1) % _Nz == 0)
				log_stream << std::endl;
		}
		log_stream << std::endl << std::endl;

		log_stream << "_zMax = " << std::endl;
		for (unsigned i = 0; i < _zMax.size(); ++i)
		{
			log_stream << _zMax[i] << '\t';

			if ((i + 1) % _Nz == 0)
				log_stream << std::endl;
		}
		log_stream << std::endl << std::endl;
	}

	void MPC_Controller::InitQP()
	{
		using namespace Eigen;

		VectorXd z_min(_Nz), z_max(_Nz);
		_platform->getAxesLimits(z_min.data(), z_max.data(), z_min.data() + _Nq, z_max.data() + _Nq, z_min.data() + 2 * _Nq, z_max.data() + 2 * _Nq);

		// Output weighting matrix
		const auto W = OutputWeightingMatrix();

		for (unsigned i = 0; i < _Nt; ++i)
		{
			// Output vector and derivatives
			VectorXd y(_Ny);
			RowMajorMatrix ssC(_Ny, _Nx);
			RowMajorMatrix ssD(_Ny, _Nu);
			_platform->Output(x(i).data(), u(i).data(), y.data(), ssC.data(), ssD.data());

			// G = [C, D]
			RowMajorMatrix G(_Ny, _Nx + _Nu);
			G << ssC, ssD;

			// H = G^T W G + \mu I
			// Adding Levenberg-Marquardt term to make H positive-definite.
			H(i) = G.transpose() * W * G + _levenbergMarquardt * MatrixXd::Identity(_Nz, _Nz);

			// g = 2 * (y_bar - y_hat) * W * G
			g(i) = 2. * (y - yRef(i)) * W * G;

			// C = [ssA, ssB];
			// x_{k+1} = C * z_k + c_k
			C(i) << getStateSpaceA(), getStateSpaceB();

			// c = K * z_k - x_{k+1}
			c(i) = C(i) * z(i) - x(i + 1);

			// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
			zMin(i) = z_min - z(i);

			// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
			zMax(i) = z_max - z(i);
		}

		H(_Nt) = _levenbergMarquardt * MatrixXd::Identity(_Nx, _Nx);
		g(_Nt).fill(0.);

		zMin(_Nt) = z_min.topRows(_Nx) - z(_Nt);
		zMax(_Nt) = z_max.topRows(_Nx) - z(_Nt);
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::H(unsigned i)
	{
		assert(i < _Nt + 1);
		const auto sz = i < _Nt ? _Nz : _Nx;
		return RowMajorMatrixMap(_H.data() + i * _Nz * _Nz, sz, sz);
	}

	MPC_Controller::VectorMap MPC_Controller::g(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_g.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::C(unsigned i)
	{
		return RowMajorMatrixMap(_C.data() + i * _Nx * _Nz, _Nx, _Nz);
	}

	MPC_Controller::VectorMap MPC_Controller::c(unsigned i)
	{
		return VectorMap(_c.data() + i * _Nx, _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::zMin(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMin.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::zMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	Eigen::MatrixXd MPC_Controller::getStateSpaceA() const
	{
		using namespace Eigen;

		MatrixXd A(_Nx, _Nx);
		A << MatrixXd::Identity(_Nq, _Nq), _sampleTime * MatrixXd::Identity(_Nq, _Nq),
			MatrixXd::Zero(_Nq, _Nq), MatrixXd::Identity(_Nq, _Nq);

		return A;
	}

	Eigen::MatrixXd MPC_Controller::getStateSpaceB() const
	{
		using namespace Eigen;

		MatrixXd B(_Nx, _Nu);
		B << _sampleTime * _sampleTime / 2. * MatrixXd::Identity(_Nq, _Nq),
			_sampleTime * MatrixXd::Identity(_Nq, _Nq);

		return B;
	}

	MPC_Controller::VectorMap MPC_Controller::z(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_z.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::x(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_z.data() + i * _Nz, _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::u(unsigned i)
	{
		assert(i < _Nt);
		return VectorMap(_z.data() + i * _Nz + _Nx, _Nu);
	}

	MPC_Controller::VectorMap MPC_Controller::xMin(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMin.data() + i * _Nz, _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::xMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::yRef(unsigned i)
	{
		assert(i < _Nt);
		return VectorMap(_yRef.data() + i * _Ny, _Ny);
	}

	void MPC_Controller::InitWorkingPoint()
	{
		using namespace Eigen;

		// Set up initial working point and reference.
		// x0 = [q0; 0]
		VectorXd x0(_Nx);
		x0.fill(0.);
		_platform->getDefaultAxesPosition(x0.data());

		// u0 = 0;
		VectorXd u0(_Nx);
		u0.fill(0.);

		VectorXd y0(_Ny);
		_platform->Output(x0.data(), u0.data(), y0.data());

		for (unsigned i = 0; i < _Nt; ++i)
			z(i) << x0, u0;

		z(_Nt) = x0;

		// Initialize QP
		InitQP();

		/** Initial MPC data setup: components not given here are set to zero (if applicable)
		*      instead of passing g, D, zLow, zUpp, one can also just pass NULL pointers (0) */
		return_t statusFlag = qpDUNES_init(&_qpData, _H.data(), _g.data(), _C.data(), _c.data(), _zMin.data(), _zMax.data(), 0, 0, 0);
		if (statusFlag != QPDUNES_OK)
			throw std::runtime_error("qpDUNES_init() failed.");
	}

}