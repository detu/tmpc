#include <MPC_Controller.hpp>

#include <stdexcept>

namespace camels
{
	MPC_Controller::MPC_Controller(unsigned state_dim, unsigned input_dim, unsigned n_path_constr, unsigned n_term_constr, double sample_time, unsigned Nt) : 
		_QP(state_dim, input_dim, n_path_constr, n_term_constr, Nt),
		_Solver(MultiStageQPSize(state_dim, input_dim, n_path_constr, n_term_constr, Nt)),
		_levenbergMarquardt(0.0),
		_sampleTime(sample_time),
		_xMin(state_dim),
		_xMax(state_dim),
		_uMin(input_dim),
		_uMax(input_dim),
		_terminalXMin(state_dim),
		_terminalXMax(state_dim),
		_Nu(input_dim),
		_Nx(state_dim),
		_Nz(input_dim + state_dim),
		_Nt(Nt),
		_Nd(n_path_constr),
		_NdT(n_term_constr)
	{
		// Allocate arrays.
		_zOpt.resize(_Nz * _Nt + _Nx);
		_w.resize(_Nz * _Nt + _Nx);

		// Initialize limits.
		_xMin.fill(-std::numeric_limits<double>::infinity());
		_xMax.fill( std::numeric_limits<double>::infinity());
		_uMin.fill(-std::numeric_limits<double>::infinity());
		_uMax.fill( std::numeric_limits<double>::infinity());
		_terminalXMin.fill(-std::numeric_limits<double>::infinity());
		_terminalXMax.fill( std::numeric_limits<double>::infinity());
	}

	void MPC_Controller::PrintQP_C(std::ostream& log_stream) const
	{
		_QP.PrintQP_C(log_stream);
	}

	void MPC_Controller::PrintQP_MATLAB(std::ostream& log_stream) const
	{
		_QP.PrintQP_MATLAB(log_stream);
	}

	void MPC_Controller::PrintWorkingPoint_MATLAB(std::ostream& os, const std::string& var_name) const
	{
		for (unsigned i = 0; i <= getNumberOfIntervals(); ++i)
			os << var_name << "{" << i + 1 << "} = [" << getWorkingPoint(i) << "];" << std::endl;
	}

	void MPC_Controller::UpdateQP()
	{
		using namespace Eigen;

		for (unsigned i = 0; i < _Nt; ++i)
			UpdateStage(i);

		MatrixXd D(_NdT, _Nx);
		VectorXd d_min(_NdT);
		VectorXd d_max(_NdT);
		TerminalConstraints(w(_Nt), D, d_min, d_max);
		_QP.D(_Nt) = D;
		_QP.dMin(_Nt) = d_min;
		_QP.dMax(_Nt) = d_max;

		_QP.zMin(_Nt) = getTerminalXMin() - w(_Nt);
		_QP.zMax(_Nt) = getTerminalXMax() - w(_Nt);
	}

	MPC_Controller::VectorMap MPC_Controller::w(unsigned i)
	{
		if(!(i < _Nt + 1))
			throw std::out_of_range("MPC_Controller::w(): index is out of range");

		return VectorMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::w(unsigned i) const
	{
		if (!(i < _Nt + 1))
			throw std::out_of_range("MPC_Controller::w(): index is out of range");

		return VectorConstMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	MPC_Controller::VectorConstMap MPC_Controller::getWorkingPoint(unsigned i) const
	{
		return w(i);
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

		// Call the QP callback, if there is one.
		if(_QPCallback)
			_QPCallback(_QP);

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

	Eigen::VectorXd MPC_Controller::getWorkingU(unsigned i) const
	{
		if (!(i < _Nt))
			throw std::out_of_range("MPC_Controller::getWorkingU(): index is out of range");

		return w(i).bottomRows(nU());
	}

	void MPC_Controller::EmbedInitialValue(const Eigen::VectorXd& x0)
	{
		// Compute linearization at new initial point.
		w(0).topRows(_Nx) = x0;
		UpdateStage(0);

		/** embed current initial value */
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

	void MPC_Controller::setXMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("MPC_Controller::setXMin(): val has a wrong size");

		_xMin = val;
	}

	void MPC_Controller::setXMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("MPC_Controller::setXMax(): val has a wrong size");

		_xMax = val;
	}

	void MPC_Controller::setTerminalXMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("MPC_Controller::setTerminalXMin(): val has a wrong size");

		_terminalXMin = val;
	}

	void MPC_Controller::setTerminalXMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("MPC_Controller::setTerminalXMax(): val has a wrong size");

		_terminalXMax = val;
	}

	void MPC_Controller::setUMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nU())
			throw std::invalid_argument("MPC_Controller::setUMin(): val has a wrong size");

		_uMin = val;
	}

	void MPC_Controller::setUMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nU())
			throw std::invalid_argument("MPC_Controller::setUMax(): val has a wrong size");

		_uMax = val;
	}

	void MPC_Controller::setQPCallback(const QPCallback& cb)
	{
		_QPCallback = cb;
	}

	void MPC_Controller::UpdateStage(unsigned i)
	{
		using namespace Eigen;

		VectorXd z_min(_Nz), z_max(_Nz);
		z_min << _xMin, _uMin;
		z_max << _xMax, _uMax;

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

		MatrixXd D(_Nd, _Nz);
		VectorXd d_min(_Nd);
		VectorXd d_max(_Nd);
		PathConstraints(i, w(i).topRows(nX()), w(i).bottomRows(nU()), D, d_min, d_max);
		_QP.D(i) = D;
		_QP.dMin(i) = d_min;
		_QP.dMax(i) = d_max;

		// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMin(i) = z_min - w(i);

		// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMax(i) = z_max - w(i);
	}
}
