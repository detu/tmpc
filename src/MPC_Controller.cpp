#include <MPC_Controller.hpp>

#include <stdexcept>

namespace rtmc
{
	MPC_Controller::MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt) 
		: _platform(platform)
		, _QP(platform->getStateDim(), platform->getInputDim(), Nt)
		, _Solver(platform->getStateDim(), platform->getInputDim(), Nt)
		, _levenbergMarquardt(0.01)
		, _sampleTime(sample_time)
		, _y(platform->getOutputDim(), Nt)
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
		_zOpt.resize(_Nz * _Nt + _Nx);
		_w.resize(_Nz * _Nt + _Nx);

		// Initialize weight matrices
		for (unsigned i = 0; i < _Nt; ++i)
			W(i).setIdentity();
	}

	MPC_Controller::~MPC_Controller()
	{
	}

	void MPC_Controller::PrintQP_C(std::ostream& log_stream) const
	{
		_QP.PrintQP_C(log_stream);
	}

	void MPC_Controller::PrintQP_MATLAB(std::ostream& log_stream) const
	{
		_QP.PrintQP_MATLAB(log_stream);
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
			_QP.H(i) = G(i).transpose() * W(i) * G(i) + _levenbergMarquardt * MatrixXd::Identity(_Nz, _Nz);

			// C = [ssA, ssB];
			// x_{k+1} = C * z_k + c_k
			RowMajorMatrix ssA(_Nx, _Nx);
			RowMajorMatrix ssB(_Nx, _Nu);
			VectorXd x_plus(_Nx);
			Integrate(w(i).data(), w(i).data() + _Nx, x_plus.data(), ssA.data(), ssB.data());
			_QP.C(i) << ssA, ssB;

			// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
			// c = f(z_k) - x_{k+1}
			_QP.c(i) = x_plus - w(i + 1).topRows(_Nx);

			// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
			_QP.zMin(i) = z_min - w(i);

			// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
			_QP.zMax(i) = z_max - w(i);
		}

		_QP.H(_Nt) = _levenbergMarquardt * MatrixXd::Identity(_Nx, _Nx);
		
		VectorXd zero_v(_Nq);
		zero_v.fill(0.);

		_QP.zMin(_Nt) << q_min, v_min;
		_QP.zMin(_Nt) -= w(_Nt);
		_QP.zMax(_Nt) << q_max, v_max;
		_QP.zMax(_Nt) -= w(_Nt);

// 		// Convert IEEE NaNs to large finite numbers to make qpDUNES happy.
// 		const double INFTY = 1.0e12;
// 
// 		for (auto& z : _zMin)
// 		{
// 			if (z == std::numeric_limits<double>::infinity())
// 				z = INFTY;
// 			if (z == -std::numeric_limits<double>::infinity())
// 				z = -INFTY;
// 		}
// 
// 		for (auto& z : _zMax)
// 		{
// 			if (z == std::numeric_limits<double>::infinity())
// 				z = INFTY;
// 			if (z == -std::numeric_limits<double>::infinity())
// 				z = -INFTY;
// 		}
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
	}

	MPC_Controller::RowMajorMatrixMap MPC_Controller::G(unsigned i)
	{
		assert(i < _Nt);
		return RowMajorMatrixMap(_G.data() + i * _Ny * _Nz, _Ny, _Nz);
	}

	void MPC_Controller::Solve()
	{
		/** solve QP */
		_Solver.Solve(_QP);
		
		// Add QP step to the working point.
		_w += _Solver.getPrimalSolution();
	}

	void MPC_Controller::PrepareForNext()
	{
		/** prepare QP for next solution */
		//qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */
		//qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

		// Shift working point
		std::copy_n(_w.data() + _Nz, (_Nt - 1) * _Nz + _Nx, _w.data());

		// Calculate new matrices.
		UpdateQP();
	}

	void MPC_Controller::getWorkingU(unsigned i, double * pu) const
	{
		assert(i < _Nt);
		std::copy_n(_w.data() + i * _Nz + _Nx, _Nu, pu);
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
			_QP.g(i) = 1. * (_y.col(i) - y_ref.col(i)).transpose() * W(i) * G(i);

		_QP.g(_Nt).setZero();
	}

	void MPC_Controller::EmbedInitialValue(const double * px0)
	{
		/** embed current initial value */
		Eigen::Map<const Eigen::VectorXd> x0(px0, _Nx);
		_QP.xMin(0) = x0 - w(0).topRows(_Nx);
		_QP.xMax(0) = x0 - w(0).topRows(_Nx);
	}
}
