#pragma once

#include <Eigen/Dense>

#include <ostream>
#include <functional>
#include <stdexcept>

namespace camels
{
	template<class _Problem, class QPSolver_>
	class ModelPredictiveController
	{
	public:
		typedef _Problem Problem;
		typedef QPSolver_ QPSolver;

		typedef typename Problem::StateVector StateVector;
		typedef typename Problem::InputVector InputVector;
		typedef typename Problem::StateInputVector StateInputVector;
		typedef typename Problem::ODEJacobianMatrix ODEJacobianMatrix;
		typedef typename Problem::LagrangeHessianMatrix LagrangeHessianMatrix;
		typedef typename Problem::MayerHessianMatrix MayerHessianMatrix;

		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;
		typedef std::function<void (typename QPSolver::MultiStageQP const&)> QPCallback;

		ModelPredictiveController(Problem const& ocp, double sample_time);
		virtual ~ModelPredictiveController();

		void InitWorkingPoint(const Eigen::VectorXd& x0);

		// Feed current state x0, get back control input u.
		InputVector Feedback(const StateVector& x0)
		{
			// Compute linearization at new initial point.
			{
				_workingPoint.w(0).topRows(_Nx) = x0;
				UpdateStage(0);

				/** embed current initial value */
				_QP.xMin(0) = x0 - _workingPoint.w(0).topRows(_Nx);
				_QP.xMax(0) = x0 - _workingPoint.w(0).topRows(_Nx);
			}

			//
			// Solve the QP.
			//
			{
				LagrangeHessianMatrix H_i;
				StateInputVector g_i;

				// Hessians and gradients of Lagrange terms.
				for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
				{
					_ocp.LagrangeTerm(i, _workingPoint.w(i), g_i, H_i);

					// Adding Levenberg-Marquardt term to make H positive-definite.
					_QP.H(i) = H_i + _levenbergMarquardt * LagrangeHessianMatrix::Identity();
					_QP.g(i) = g_i;
				}

				// Hessian and gradient of Mayer term.
				typename Problem::MayerHessianMatrix H_T;
				typename Problem::StateVector g_T;
				_ocp.MayerTerm(_workingPoint.wend(), g_T, H_T);

				// Adding Levenberg-Marquardt term to make H positive-definite.
				_QP.Hend() = H_T + _levenbergMarquardt * MayerHessianMatrix::Identity();
				_QP.gend() = g_T;

				// Call the QP callback, if there is one.
				if(_QPCallback)
					_QPCallback(_QP);

				/** solve QP */
				_Solver.Solve(_QP, _solution);
			}

			// Return the calculated control input.
			return _workingPoint.w(0).bottomRows(nU()) + _solution.w(0).template bottomRows<_Nu>();
		}

		void Preparation()
		{
			// Add QP step to the working point.
			_workingPoint += _solution;

			/** prepare QP for next solution */
			//qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */
			//qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

			// Shift working point
			_workingPoint.shift();

			// Calculate new matrices.
			UpdateQP();
		}

		void PrintQP_C(std::ostream& os) const;
		void PrintQP_zMax_C(std::ostream& log_stream) const;
		void PrintQP_zMin_C(std::ostream& log_stream) const;
		void PrintQP_MATLAB(std::ostream& log_stream) const;

		// Log working point
		void PrintWorkingPoint_MATLAB(std::ostream& os, const std::string& var_name) const;
		
		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		double getSampleTime() const;

		unsigned nU() const;
		unsigned nX() const;
		unsigned nZ() const { return _Nz; }

		void setQPCallback(const QPCallback& cb);

	private:

		// Initializes _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		void UpdateStage(unsigned i);

		// Private data members.
		Problem const& _ocp;

		const double _sampleTime;

		static const unsigned _Nu = Problem::NU;
		static const unsigned _Nx = Problem::NX;
		static const unsigned _Nz = Problem::NW;
		static const unsigned _Nd = Problem::NC;
		static const unsigned _NdT = Problem::NCT;
		
		typename QPSolver::MultiStageQP _QP;
		typename QPSolver::Point _solution;
		QPSolver _Solver;

		// A callback to call back before solving each QP.
		QPCallback _QPCallback;

		double _levenbergMarquardt;

		// Working point (linearization point).
		typename QPSolver::Point _workingPoint;
	};

	template<class _Problem, class QPSolver_>
	ModelPredictiveController<_Problem, QPSolver_>::ModelPredictiveController(Problem const& ocp, double sample_time)
	:	_ocp(ocp)
	,	_QP(ocp.getNumberOfIntervals())
	,	_workingPoint(ocp.getNumberOfIntervals())
	,	_solution(ocp.getNumberOfIntervals())
	,	_Solver(ocp.getNumberOfIntervals()),
		_levenbergMarquardt(0.0),
		_sampleTime(sample_time)
	{
	}

	template<class _Problem, class QPSolver_>
	ModelPredictiveController<_Problem, QPSolver_>::~ModelPredictiveController()
	{
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::PrintQP_C(std::ostream& log_stream) const
	{
		_QP.PrintQP_C(log_stream);
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::PrintQP_MATLAB(std::ostream& log_stream) const
	{
		_QP.PrintQP_MATLAB(log_stream);
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::PrintWorkingPoint_MATLAB(std::ostream& os, const std::string& var_name) const
	{
		for (unsigned i = 0; i < nT(); ++i)
			os << var_name << "{" << i + 1 << "} = [" << _workingPoint.w(i) << "];" << std::endl;
		os << var_name << "{" << nT() + 1 << "} = [" << _workingPoint.wend() << "];" << std::endl;
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::UpdateQP()
	{
		for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
			UpdateStage(i);

		typename Problem::TerminalConstraintJacobianMatrix D;
		typename Problem::TerminalConstraintVector d_min, d_max;
		_ocp.TerminalConstraints(_workingPoint.wend(), D, d_min, d_max);
		_QP.Dend() = D;
		_QP.dendMin() = d_min;
		_QP.dendMax() = d_max;

		_QP.zendMin() = _ocp.getTerminalStateMin() - _workingPoint.wend();
		_QP.zendMax() = _ocp.getTerminalStateMax() - _workingPoint.wend();
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::InitWorkingPoint( const Eigen::VectorXd& x0 )
	{
		// u0 = 0;
		InputVector const u0 = InputVector::Zero();

		for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
			_workingPoint.w(i) << x0, u0;

		_workingPoint.wend() = x0;

		// Initialize QP
		UpdateQP();
	}

	template<class _Problem, class QPSolver_>
	double ModelPredictiveController<_Problem, QPSolver_>::getSampleTime() const
	{
		return _sampleTime;
	}

	template<class _Problem, class QPSolver_>
	unsigned ModelPredictiveController<_Problem, QPSolver_>::nU() const
	{
		return _Nu;
	}

	template<class _Problem, class QPSolver_>
	unsigned ModelPredictiveController<_Problem, QPSolver_>::nX() const
	{
		return _Nx;
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::setQPCallback(const QPCallback& cb)
	{
		_QPCallback = cb;
	}

	template<class _Problem, class QPSolver_>
	void ModelPredictiveController<_Problem, QPSolver_>::UpdateStage(unsigned i)
	{
		StateInputVector z_min, z_max;
		z_min << _ocp.getStateMin(), _ocp.getInputMin();
		z_max << _ocp.getStateMax(), _ocp.getInputMax();

		// C = [ssA, ssB];
		// x_{k+1} = C * z_k + c_k
		typename Problem::StateVector x_plus;
		typename Problem::ODEJacobianMatrix J;
		_ocp.Integrate(_workingPoint.w(i), getSampleTime(), x_plus, J);
		_QP.C(i) = J;

		// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
		// c = f(z_k) - x_{k+1}
		_QP.c(i) = x_plus - _workingPoint.x(i + 1);

		typename Problem::ConstraintJacobianMatrix D;
		typename Problem::ConstraintVector d_min, d_max;
		_ocp.PathConstraints(i, _workingPoint.w(i), D, d_min, d_max);
		_QP.D(i) = D;
		_QP.dMin(i) = d_min;
		_QP.dMax(i) = d_max;

		// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMin(i) = z_min - _workingPoint.w(i);

		// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMax(i) = z_max - _workingPoint.w(i);
	}
}
