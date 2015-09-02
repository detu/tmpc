#include <MPC_Controller.hpp>

#include <stdexcept>

namespace camels
{
	MPC_Controller::MPC_Controller(unsigned state_dim, unsigned input_dim, double sample_time, unsigned Nt) 
		: _QP(state_dim, input_dim, Nt)
		, _Solver(state_dim, input_dim, Nt)
		, _levenbergMarquardt(0.0)
		, _sampleTime(sample_time)
		, _xMin(state_dim)
		, _xMax(state_dim)
		, _uMin(input_dim)
		, _uMax(input_dim)
	{
		// Get sizes.
		_Nu = input_dim;
		_Nx = state_dim;
		_Nz = _Nx + _Nu;
		_Nt = Nt;

		// Allocate arrays.
		_zOpt.resize(_Nz * _Nt + _Nx);
		_w.resize(_Nz * _Nt + _Nx);

		// Initialize limits.
		_xMin.fill(-std::numeric_limits<double>::infinity());
		_xMax.fill( std::numeric_limits<double>::infinity());
		_uMin.fill(-std::numeric_limits<double>::infinity());
		_uMax.fill( std::numeric_limits<double>::infinity());
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
		z_min << _xMin, _uMin;
		z_max << _xMax, _uMax;
		
		/*
		const auto q_min = z_min.middleRows(      0, _Nq), q_max = z_max.middleRows(      0, _Nq);
		const auto v_min = z_min.middleRows(    _Nq, _Nq), v_max = z_max.middleRows(    _Nq, _Nq);
		const auto u_min = z_min.middleRows(2 * _Nq, _Nq), u_max = z_max.middleRows(2 * _Nq, _Nq);
		const VectorXd q_min_final = q_min + v_min.cwiseAbs2().cwiseQuotient(2 * u_max);
		const VectorXd q_max_final = q_max + v_max.cwiseAbs2().cwiseQuotient(2 * u_min);
		*/

		for (unsigned i = 0; i < _Nt; ++i)
		{
			// C = [ssA, ssB];
			// x_{k+1} = C * z_k + c_k
			MatrixXd ssA(_Nx, _Nx);
			MatrixXd ssB(_Nx, _Nu);
			VectorXd x_plus(_Nx);
			Integrate(w(i).topRows(nX()), w(i).bottomRows(nU()), x_plus, ssA, ssB);
			_QP.C(i) << ssA, ssB;

			// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
			// c = f(z_k) - x_{k+1}
			_QP.c(i) = x_plus - w(i + 1).topRows(_Nx);

			// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
			_QP.zMin(i) = z_min - w(i);

			// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
			_QP.zMax(i) = z_max - w(i);
		}

		_QP.zMin(_Nt) = getXMin() - w(_Nt);
		_QP.zMax(_Nt) = getXMax() - w(_Nt);
	}

	MPC_Controller::VectorMap MPC_Controller::w(unsigned i)
	{
		assert(i < _Nt + 1);
		return VectorMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	void MPC_Controller::InitWorkingPoint( const Eigen::VectorXd& x0 )
	{
		using namespace Eigen;

		// u0 = 0;
		VectorXd u0(_Nu);
		u0.fill(0.);

		for (unsigned i = 0; i < _Nt; ++i)
			w(i) << x0, u0;

		w(_Nt) = x0;

		// Initialize QP
		UpdateQP();
	}

	void MPC_Controller::Solve()
	{
		using Eigen::MatrixXd;
		using Eigen::VectorXd;

		MatrixXd H_i(nZ(), nZ());
		VectorXd g_i(nZ());

		// Hessians and gradients of Lagrange terms.
		for (unsigned i = 0; i < _Nt; ++i)
		{
			LagrangeTerm(w(i), i, H_i, g_i);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			_QP.H(i) = H_i + _levenbergMarquardt * MatrixXd::Identity(_Nz, _Nz);
			_QP.g(i) = g_i;
		}

		// Hessian and gradient of Mayer term.
		MatrixXd H_T(nX(), nX());
		VectorXd g_T(nX());
		MayerTerm(w(_Nt), H_T, g_T);

		// Adding Levenberg-Marquardt term to make H positive-definite.
		_QP.H(_Nt) = H_T + _levenbergMarquardt * MatrixXd::Identity(_Nx, _Nx);
		_QP.g(_Nt) = g_T;

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

	void MPC_Controller::EmbedInitialValue(const double * px0)
	{
		/** embed current initial value */
		Eigen::Map<const Eigen::VectorXd> x0(px0, _Nx);
		_QP.xMin(0) = x0 - w(0).topRows(_Nx);
		_QP.xMax(0) = x0 - w(0).topRows(_Nx);
	}

	double MPC_Controller::getSampleTime() const
	{
		return _sampleTime;
	}

	unsigned MPC_Controller::nU() const
	{
		return _Nu;
	}

	unsigned MPC_Controller::nX() const
	{
		return _Nx;
	}
}
