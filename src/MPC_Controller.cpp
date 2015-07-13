#include <MPC_Controller.hpp>

#include <stdexcept>

std::string qpDUNES_getErrorMessage(return_t code);

class qpDUNESException : public std::runtime_error
{
public:
	qpDUNESException(return_t code, const std::string& info)
		: std::runtime_error(qpDUNES_getErrorMessage(code))
		, _code(code)
	{
		
	}

private:
	return_t _code;
};

namespace rtmc
{
	MPC_Controller::MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, 
		unsigned Nt, const qpOptions_t * qpOptions /*= nullptr*/, std::ostream * call_log /*= nullptr*/) : _platform(platform)
		, _levenbergMarquardt(0.01)
		, _sampleTime(sample_time)
		, _y(platform->getOutputDim(), Nt)
		, _qpDUNES_call_log(call_log)
	{
		// Get sizes.
		_Nq = platform->getNumberOfAxes();
		_Nu = platform->getInputDim();
		_Nx = platform->getStateDim();
		_Ny = platform->getOutputDim();
		_Nz = _Nx + _Nu;
		_Nt = Nt;

		// Allocate arrays.
		_W.resize(_Ny * _Ny * _Nt);
		_G.resize(_Ny * _Nz * _Nt);
		_H.resize(_Nz * _Nz * _Nt + _Nx * _Nx);
		_g.resize(_Nz * _Nt + _Nx);
		_C.resize(_Nx * _Nz * _Nt);
		_c.resize(_Nx * _Nt);
		_zMin.resize(_Nz * _Nt + _Nx);
		_zMax.resize(_Nz * _Nt + _Nx);
		_zOpt.resize(_Nz * _Nt + _Nx);
		_w.resize(_Nz * _Nt + _Nx);

		// Initialize weight matrices
		for (unsigned i = 0; i < _Nt; ++i)
			W(i).setIdentity();

		// Set up qpData
		unsigned int* nD = 0;	  			/* number of affine constraints */
		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			*_qpDUNES_call_log << "const unsigned int nT = " << _Nt << ";" << std::endl;
			*_qpDUNES_call_log << "const unsigned int nX = " << _Nx << ";" << std::endl;
			*_qpDUNES_call_log << "const unsigned int nU = " << _Nu << ";" << std::endl;
			*_qpDUNES_call_log << "unsigned int * nD = " << nD << ";" << std::endl;
			*_qpDUNES_call_log << "qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();" << std::endl;
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_setup(&qpData, nT, nX, nU, nD, &qpOptions);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}

		return_t statusFlag = qpDUNES_setup(&_qpData, _Nt, _Nx, _Nu, nD, qpOptions);
		if (statusFlag != QPDUNES_OK)
			throw qpDUNESException(statusFlag, "qpDUNES_setup() failed.");
	}

	MPC_Controller::~MPC_Controller()
	{
		/** cleanup of allocated data */
		if (_qpDUNES_call_log)
			*_qpDUNES_call_log << "qpDUNES_cleanup(&qpData);";

		qpDUNES_cleanup(&_qpData);
	}

	void MPC_Controller::PrintQP_C(std::ostream& log_stream) const
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

	void MPC_Controller::PrintQP_MATLAB(std::ostream& log_stream) const
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

	void MPC_Controller::UpdateQP()
	{
		using namespace Eigen;

		VectorXd z_min(_Nz), z_max(_Nz);
		_platform->getAxesLimits(z_min.data(), z_max.data(), z_min.data() + _Nq, z_max.data() + _Nq, z_min.data() + 2 * _Nq, z_max.data() + 2 * _Nq);
		const auto q_min = z_min.middleRows(      0, _Nq), q_max = z_max.middleRows(      0, _Nq);
		const auto v_min = z_min.middleRows(    _Nq, _Nq), v_max = z_max.middleRows(    _Nq, _Nq);
		const auto u_min = z_min.middleRows(2 * _Nq, _Nq), u_max = z_max.middleRows(2 * _Nq, _Nq);
		const VectorXd q_min_final = q_min + v_min.cwiseAbs2().cwiseQuotient(2 * u_max);
		const VectorXd q_max_final = q_max + v_max.cwiseAbs2().cwiseQuotient(2 * u_min);

		for (unsigned i = 0; i < _Nt; ++i)
		{
			// Output vector and derivatives.
			// Output() returns derivative matrices in column-major format.
			MatrixXd ssC(_Ny, _Nx);
			MatrixXd ssD(_Ny, _Nu);
			_platform->Output(w(i).data(), w(i).data() + _Nx, _y.col(i).data(), ssC.data(), ssD.data());

			// G = [C, D]
			G(i) << ssC, ssD;

			// H = G^T W G + \mu I
			// Adding Levenberg-Marquardt term to make H positive-definite.
			H(i) = G(i).transpose() * W(i) * G(i) + _levenbergMarquardt * MatrixXd::Identity(_Nz, _Nz);

			// C = [ssA, ssB];
			// x_{k+1} = C * z_k + c_k
			RowMajorMatrix ssA(_Nx, _Nx);
			RowMajorMatrix ssB(_Nx, _Nu);
			VectorXd x_plus(_Nx);
			Integrate(w(i).data(), w(i).data() + _Nx, x_plus.data(), ssA.data(), ssB.data());
			C(i) << ssA, ssB;

			// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
			// c = f(z_k) - x_{k+1}
			c(i) = x_plus - w(i + 1).topRows(_Nx);

			// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
			zMin(i) = z_min - w(i);

			// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
			zMax(i) = z_max - w(i);
		}

		H(_Nt) = _levenbergMarquardt * MatrixXd::Identity(_Nx, _Nx);
		
		VectorXd zero_v(_Nq);
		zero_v.fill(0.);

		zMin(_Nt) << q_min, v_min;
		zMin(_Nt) -= w(_Nt);
		zMax(_Nt) << q_max, v_max;
		zMax(_Nt) -= w(_Nt);

		// Convert IEEE NaNs to large finite numbers to make qpDUNES happy.
		const double INFTY = 1.0e12;
		
		for (auto& z : _zMin)
		{
			if (z == std::numeric_limits<double>::infinity())
				z = INFTY;
			if (z == -std::numeric_limits<double>::infinity())
				z = -INFTY;
		}

		for (auto& z : _zMax)
		{
			if (z == std::numeric_limits<double>::infinity())
				z = INFTY;
			if (z == -std::numeric_limits<double>::infinity())
				z = -INFTY;
		}
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::H(unsigned i)
	{
		assert(i < _Nt + 1);
		const auto sz = i < _Nt ? _Nz : _Nx;
		return RowMajorMatrixMap(_H.data() + i * _Nz * _Nz, sz, sz);
	}

	MPC_Controller::RowMajorMatrixConstMap MPC_Controller::H(unsigned i) const
	{
		assert(i < _Nt + 1);
		const auto sz = i < _Nt ? _Nz : _Nx;
		return RowMajorMatrixConstMap(_H.data() + i * _Nz * _Nz, sz, sz);
	}

	MPC_Controller::VectorMap MPC_Controller::g(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_g.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::g(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_g.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::C(unsigned i)
	{
		return RowMajorMatrixMap(_C.data() + i * _Nx * _Nz, _Nx, _Nz);
	}

	MPC_Controller::RowMajorMatrixConstMap MPC_Controller::C(unsigned i) const
	{
		return RowMajorMatrixConstMap(_C.data() + i * _Nx * _Nz, _Nx, _Nz);
	}

	MPC_Controller::VectorMap MPC_Controller::c(unsigned i)
	{
		return VectorMap(_c.data() + i * _Nx, _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::c(unsigned i) const
	{
		return VectorConstMap(_c.data() + i * _Nx, _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::zMin(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMin.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::zMin(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMin.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorMap MPC_Controller::zMax(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_zMax.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::zMax(unsigned i) const
	{
		assert(i < _Nt + 1);
		return VectorConstMap(_zMax.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
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

	MPC_Controller::VectorMap MPC_Controller::w(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
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

	void MPC_Controller::InitWorkingPoint()
	{
		using namespace Eigen;

		// Set up initial working point and reference.
		// x0 = [q0; 0]
		VectorXd x0(_Nx);
		x0.fill(0.);
		_platform->getDefaultAxesPosition(x0.data());

		// u0 = 0;
		VectorXd u0(_Nu);
		u0.fill(0.);

		VectorXd y0(_Ny);
		_platform->Output(x0.data(), u0.data(), y0.data());

		for (unsigned i = 0; i < _Nt; ++i)
			w(i) << x0, u0;

		w(_Nt) = x0;

		// Initialize QP
		UpdateQP();

// 		/** set sparsity of primal Hessian and local constraint matrix */
// 		for (unsigned kk = 0; kk < _Nt + 1; ++kk) {
// 			_qpData.intervals[kk]->H.sparsityType = QPDUNES_DENSE;
// 			//_qpData.intervals[kk]->D.sparsityType = QPDUNES_IDENTITY;
// 		}

		/** Initial MPC data setup: components not given here are set to zero (if applicable)
		*      instead of passing g, D, zLow, zUpp, one can also just pass NULL pointers (0) */
		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			PrintQP_C(*_qpDUNES_call_log);
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_init(&qpData, H, g, C, c, zLow, zUpp, 0, 0, 0);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}

		return_t statusFlag = qpDUNES_init(&_qpData, _H.data(), _g.data(), _C.data(), _c.data(), _zMin.data(), _zMax.data(), 0, 0, 0);
		if (statusFlag != QPDUNES_OK)
			throw qpDUNESException(statusFlag, "qpDUNES_init() failed.");
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::G(unsigned i)
	{
		assert(i < _Nt);
		return RowMajorMatrixMap(_G.data() + i * _Ny * _Nz, _Ny, _Nz);
	}

	void MPC_Controller::Solve()
{
		/** solve QP */
		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_solve(&qpData);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}

		return_t statusFlag = qpDUNES_solve(&_qpData);
		if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND)
			throw qpDUNESException(statusFlag, "QP solution failed (qpDUNES_solve()).");

		/** obtain optimal solution */
		qpDUNES_getPrimalSol(&_qpData, _zOpt.data());
		//qpDUNES_getDualSol(&qpData, lambdaOpt, muOpt);
		
		// Add QP step to the working point.
		std::transform(_w.begin(), _w.end(), _zOpt.begin(), _w.begin(), std::plus<double>());
	}

	void MPC_Controller::PrepareForNext()
	{
		/** prepare QP for next solution */
		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_shiftLambda(&qpData);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}
		qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */

		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_shiftIntervals(&qpData);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}
		qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

		// Shift working point
		std::copy_n(_w.begin() + _Nz, (_Nt - 1) * _Nz + _Nx, _w.begin());

		// Calculate new matrices.
		UpdateQP();
	}

	void MPC_Controller::getWorkingU(unsigned i, double * pu) const
	{
		assert(i < _Nt);
		std::copy_n(_w.begin() + i * _Nz + _Nx, _Nu, pu);
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::W(unsigned i)
	{
		assert(i < _Nt);
		return RowMajorMatrixMap(_W.data() + i * _Ny * _Ny, _Ny, _Ny);
	}

	void MPC_Controller::Integrate(const double * px, const double * pu, double * px_next, double * pA, double * pB) const
	{
		using namespace Eigen;

		Map<const VectorXd> x(px, _Nx);
		Map<const VectorXd> u(pu, _Nu);
		Map<VectorXd> x_next(px_next, _Nx);
		RowMajorMatrixMap A(pA, _Nx, _Nx);
		RowMajorMatrixMap B(pB, _Nx, _Nu);

		A = getStateSpaceA();
		B = getStateSpaceB();

		x_next = A * x + B * u;
	}

	void MPC_Controller::SetReference(const double * py_ref)
	{
		// Update g
		Eigen::Map<const Eigen::MatrixXd> y_ref(py_ref, _Ny, _Nt);

		for (unsigned i = 0; i < _Nt; ++i)
			// g = 2 * (y_bar - y_hat)^T * W * G
			g(i) = 1. * (_y.col(i) - y_ref.col(i)).transpose() * W(i) * G(i);

		g(_Nt).fill(0.);

		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			PrintQP_C(*_qpDUNES_call_log);
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_updateData(&qpData, H, g, C, c, zLow, zUpp, 0, 0, 0);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}

		return_t statusFlag = qpDUNES_updateData(&_qpData, _H.data(), _g.data(), _C.data(), _c.data(), _zMin.data(), _zMax.data(), 0, 0, 0);		/* data update: components not given here keep their previous value */
		if (statusFlag != QPDUNES_OK)
			throw qpDUNESException(statusFlag, "Data update failed (qpDUNES_updateData()).");
	}

	void MPC_Controller::EmbedInitialValue(const double * px0)
	{
		/** embed current initial value */
		Eigen::Map<const Eigen::VectorXd> x0(px0, _Nx);
		zMin(0).topRows(_Nx) = x0 - w(0).topRows(_Nx);
		zMax(0).topRows(_Nx) = x0 - w(0).topRows(_Nx);

		if (_qpDUNES_call_log)
		{
			*_qpDUNES_call_log << "{" << std::endl;
			PrintQP_zMin_C(*_qpDUNES_call_log);
			PrintQP_zMax_C(*_qpDUNES_call_log);
			*_qpDUNES_call_log << "return_t statusFlag = qpDUNES_updateIntervalData(&qpData, qpData.intervals[0], 0, 0, 0, 0, zLow, zUpp, 0, 0, 0, 0);" << std::endl;
			*_qpDUNES_call_log << "assert(statusFlag == QPDUNES_OK);" << std::endl;
			*_qpDUNES_call_log << "}" << std::endl;
		}

		return_t statusFlag = qpDUNES_updateIntervalData(&_qpData, _qpData.intervals[0], 0, 0, 0, 0, zMin(0).data(), zMax(0).data(), 0, 0, 0, 0);
		if (statusFlag != QPDUNES_OK)
			throw qpDUNESException(statusFlag, "Initial value embedding failed (qpDUNES_updateIntervalData()).");
	}

	void MPC_Controller::PrintQP_zMin_C(std::ostream& log_stream) const
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

	void MPC_Controller::PrintQP_zMax_C(std::ostream& log_stream) const
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

#define QPDUNES_ERROR_MESSAGE_CASE( k ) case k: return (#k);

std::string qpDUNES_getErrorMessage(return_t code)
{
	switch (code)
	{
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_OK)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_SUCC_SUBOPTIMAL_TERMINATION)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_STAGE_QP_INFEASIBLE)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_STAGE_COUPLING_INFEASIBLE)

		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_UNKNOWN_ERROR)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_UNKNOWN_MATRIX_SPARSITY_TYPE)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_UNKNOWN_LS_TYPE)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_INVALID_ARGUMENT)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_ITERATION_LIMIT_REACHED)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_DIVISION_BY_ZERO)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_NUMBER_OF_MAX_LINESEARCH_ITERATIONS_REACHED)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_DECEEDED_MIN_LINESEARCH_STEPSIZE)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_EXCEEDED_MAX_LINESEARCH_STEPSIZE)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_ERR_NEWTON_SYSTEM_NO_ASCENT_DIRECTION)
		QPDUNES_ERROR_MESSAGE_CASE(QPDUNES_NOTICE_NEWTON_MATRIX_NOT_SET_UP)

	default:
		return "Unknown qpDUNES error code " + std::to_string(code);
	}
}