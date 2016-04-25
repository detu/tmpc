#pragma once

#include <MultiStageQP.hpp>
#include <CondensingSolver.hpp>

#include <Eigen/Dense>

#include <vector>
#include <memory>
#include <ostream>
#include <functional>

namespace camels
{
	template<class _Problem>
	class ModelPredictiveController
	{
	public:
		typedef _Problem Problem;
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
		typedef std::function<void (const MultiStageQP&)> QPCallback;

		ModelPredictiveController(Problem const& ocp, double sample_time);
		virtual ~ModelPredictiveController();

		void InitWorkingPoint(const Eigen::VectorXd& x0);
		void Solve();

		void EmbedInitialValue(const Eigen::VectorXd& x0);

		Eigen::VectorXd getWorkingU(unsigned i) const;
		void PrepareForNext();

		void PrintQP_C(std::ostream& os) const;
		void PrintQP_zMax_C(std::ostream& log_stream) const;
		void PrintQP_zMin_C(std::ostream& log_stream) const;
		void PrintQP_MATLAB(std::ostream& log_stream) const;

		// Log working point
		void PrintWorkingPoint_MATLAB(std::ostream& os, const std::string& var_name) const;
		
		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		unsigned getNumberOfIntervals() const { return _ocp.getNumberOfIntervals(); }
		double getSampleTime() const;

		unsigned nU() const;
		unsigned nX() const;
		unsigned nZ() const { return _Nz; }

		void setQPCallback(const QPCallback& cb);

		VectorConstMap getWorkingPoint(unsigned i) const;

	protected:
		// TODO: get rid of these functions, call directly from _ocp.
		const Eigen::VectorXd getXMin() const { return _ocp.getStateMin(); }
		const Eigen::VectorXd getXMax() const { return _ocp.getStateMax(); }
		const Eigen::VectorXd getTerminalXMin() const { return _ocp.getStateMin(); }
		const Eigen::VectorXd getTerminalXMax() const { return _ocp.getStateMax(); }
		const Eigen::VectorXd getUMin() const { return _ocp.getInputMin(); }
		const Eigen::VectorXd getUMax() const { return _ocp.getInputMax(); }

	private:
		// Initialized _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		void UpdateStage(unsigned i);

		VectorMap w(unsigned i);
		VectorConstMap w(unsigned i) const;

		// Private data membets.
		Problem const& _ocp;

		const double _sampleTime;

		static const unsigned _Nu = Problem::NU;
		static const unsigned _Nx = Problem::NX;
		static const unsigned _Nz = Problem::NW;
		static const unsigned _Nd = Problem::NC;
		static const unsigned _NdT = Problem::NCT;
		
		MultiStageQP _QP;
		CondensingSolver _Solver;

		// A callback to call back before solving each QP.
		QPCallback _QPCallback;

		double _levenbergMarquardt;

		// Primal optimal solution.
		// _zOpt stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zOpt;

		// Working point (linearization point).
		// _w stores _Nt vectors of size _Nz and 1 vector of size _Nx
		Eigen::VectorXd _w;

	};

	template<class _Problem>
	ModelPredictiveController<_Problem>::ModelPredictiveController(Problem const& ocp, double sample_time)
	:	_ocp(ocp)
	,	_QP(_Nx, _Nu, _Nd, _NdT, ocp.getNumberOfIntervals()),
		_Solver(MultiStageQPSize(_Nx, _Nu, _Nd, _NdT, ocp.getNumberOfIntervals())),
		_levenbergMarquardt(0.0),
		_sampleTime(sample_time)
	{
		// Allocate arrays.
		_zOpt.resize(_Nz * ocp.getNumberOfIntervals() + _Nx);
		_w.resize(_Nz * ocp.getNumberOfIntervals() + _Nx);
	}

	template<class _Problem>
	ModelPredictiveController<_Problem>::~ModelPredictiveController()
	{
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::PrintQP_C(std::ostream& log_stream) const
	{
		_QP.PrintQP_C(log_stream);
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::PrintQP_MATLAB(std::ostream& log_stream) const
	{
		_QP.PrintQP_MATLAB(log_stream);
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::PrintWorkingPoint_MATLAB(std::ostream& os, const std::string& var_name) const
	{
		for (unsigned i = 0; i <= getNumberOfIntervals(); ++i)
			os << var_name << "{" << i + 1 << "} = [" << getWorkingPoint(i) << "];" << std::endl;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::UpdateQP()
	{
		using namespace Eigen;

		for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
			UpdateStage(i);

		typename Problem::TerminalConstraintJacobianMatrix D;
		typename Problem::TerminalConstraintVector d_min, d_max;
		_ocp.TerminalConstraints(w(_ocp.getNumberOfIntervals()), D, d_min, d_max);
		_QP.D(getNumberOfIntervals()) = D;
		_QP.dMin(getNumberOfIntervals()) = d_min;
		_QP.dMax(getNumberOfIntervals()) = d_max;

		_QP.zMin(getNumberOfIntervals()) = getTerminalXMin() - w(getNumberOfIntervals());
		_QP.zMax(getNumberOfIntervals()) = getTerminalXMax() - w(getNumberOfIntervals());
	}

	template<class _Problem>
	typename ModelPredictiveController<_Problem>::VectorMap ModelPredictiveController<_Problem>::w(unsigned i)
	{
		if(!(i < _ocp.getNumberOfIntervals() + 1))
			throw std::out_of_range("ModelPredictiveController<_Problem>::w(): index is out of range");

		return VectorMap(_w.data() + i * _Nz, i < _ocp.getNumberOfIntervals() ? _Nz : _Nx);
	}

	template<class _Problem>
	typename ModelPredictiveController<_Problem>::VectorConstMap ModelPredictiveController<_Problem>::w(unsigned i) const
	{
		if (!(i < _ocp.getNumberOfIntervals() + 1))
			throw std::out_of_range("ModelPredictiveController<_Problem>::w(): index is out of range");

		return VectorConstMap(_w.data() + i * _Nz, i < _ocp.getNumberOfIntervals() ? _Nz : _Nx);
	}

	template<class _Problem>
	typename ModelPredictiveController<_Problem>::VectorConstMap ModelPredictiveController<_Problem>::getWorkingPoint(unsigned i) const
	{
		return w(i);
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::InitWorkingPoint( const Eigen::VectorXd& x0 )
	{
		using namespace Eigen;

		// u0 = 0;
		InputVector const u0 = InputVector::Zero();

		for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
			w(i) << x0, u0;

		w(_ocp.getNumberOfIntervals()) = x0;

		// Initialize QP
		UpdateQP();
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::Solve()
	{
		using Eigen::MatrixXd;
		using Eigen::VectorXd;

		LagrangeHessianMatrix H_i;
		StateInputVector g_i;

		// Hessians and gradients of Lagrange terms.
		for (unsigned i = 0; i < _ocp.getNumberOfIntervals(); ++i)
		{
			_ocp.LagrangeTerm(i, w(i), g_i, H_i);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			// TODO: Move it to the OCP class.
			_QP.H(i) = H_i + _levenbergMarquardt * MatrixXd::Identity(_Nz, _Nz);
			_QP.g(i) = g_i;
		}

		// Hessian and gradient of Mayer term.
		typename Problem::MayerHessianMatrix H_T;
		typename Problem::StateVector g_T;
		_ocp.MayerTerm(w(_ocp.getNumberOfIntervals()), g_T, H_T);

		// Adding Levenberg-Marquardt term to make H positive-definite.
		_QP.H(getNumberOfIntervals()) = H_T + _levenbergMarquardt * MayerHessianMatrix::Identity();
		_QP.g(getNumberOfIntervals()) = g_T;

		// Call the QP callback, if there is one.
		if(_QPCallback)
			_QPCallback(_QP);

		/** solve QP */
		_Solver.Solve(_QP);

		// Add QP step to the working point.
		_w += _Solver.getPrimalSolution();
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::PrepareForNext()
	{
		/** prepare QP for next solution */
		//qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */
		//qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

		// Shift working point
		std::copy_n(_w.data() + _Nz, (_ocp.getNumberOfIntervals() - 1) * _Nz + _Nx, _w.data());

		// Calculate new matrices.
		UpdateQP();
	}

	template<class _Problem>
	Eigen::VectorXd ModelPredictiveController<_Problem>::getWorkingU(unsigned i) const
	{
		if (!(i < _ocp.getNumberOfIntervals()))
			throw std::out_of_range("ModelPredictiveController<_Problem>::getWorkingU(): index is out of range");

		return w(i).bottomRows(nU());
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::EmbedInitialValue(const Eigen::VectorXd& x0)
	{
		// Compute linearization at new initial point.
		w(0).topRows(_Nx) = x0;
		UpdateStage(0);

		/** embed current initial value */
		_QP.xMin(0) = x0 - w(0).topRows(_Nx);
		_QP.xMax(0) = x0 - w(0).topRows(_Nx);
	}

	template<class _Problem>
	double ModelPredictiveController<_Problem>::getSampleTime() const
	{
		return _sampleTime;
	}

	template<class _Problem>
	unsigned ModelPredictiveController<_Problem>::nU() const
	{
		return _Nu;
	}

	template<class _Problem>
	unsigned ModelPredictiveController<_Problem>::nX() const
	{
		return _Nx;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setQPCallback(const QPCallback& cb)
	{
		_QPCallback = cb;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::UpdateStage(unsigned i)
	{
		StateInputVector z_min, z_max;
		z_min << _ocp.getStateMin(), _ocp.getInputMin();
		z_max << _ocp.getStateMax(), _ocp.getInputMax();

		// C = [ssA, ssB];
		// x_{k+1} = C * z_k + c_k
		typename Problem::StateVector x_plus;
		typename Problem::ODEJacobianMatrix J;
		_ocp.Integrate(w(i), getSampleTime(), x_plus, J);
		_QP.C(i) = J;

		// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
		// c = f(z_k) - x_{k+1}
		_QP.c(i) = x_plus - w(i + 1).topRows(_Nx);

		typename Problem::ConstraintJacobianMatrix D;
		typename Problem::ConstraintVector d_min, d_max;
		_ocp.PathConstraints(i, w(i), D, d_min, d_max);
		_QP.D(i) = D;
		_QP.dMin(i) = d_min;
		_QP.dMax(i) = d_max;

		// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMin(i) = z_min - w(i);

		// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
		_QP.zMax(i) = z_max - w(i);
	}
}
