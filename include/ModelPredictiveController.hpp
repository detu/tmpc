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
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;
		typedef std::function<void (const MultiStageQP&)> QPCallback;

		ModelPredictiveController(unsigned state_dim, unsigned input_dim, unsigned n_path_constr, unsigned n_term_constr, double sample_time, unsigned Nt);
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

		unsigned getNumberOfIntervals() const { return _Nt; }
		double getSampleTime() const;

		unsigned nU() const;
		unsigned nX() const;
		unsigned nZ() const { return _Nz; }

		const Eigen::VectorXd& getXMin() const { return _xMin; }
		void setXMin(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getXMax() const { return _xMax; }
		void setXMax(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getTerminalXMin() const { return _terminalXMin; }
		void setTerminalXMin(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getTerminalXMax() const { return _terminalXMax; }
		void setTerminalXMax(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getUMin() const { return _uMin; }
		void setUMin(const Eigen::VectorXd& val);

		const Eigen::VectorXd& getUMax() const { return _uMax; }
		void setUMax(const Eigen::VectorXd& val);

		void setQPCallback(const QPCallback& cb);

		VectorConstMap getWorkingPoint(unsigned i) const;

	protected:
		virtual void LagrangeTerm(const Eigen::VectorXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) const = 0;
		virtual void MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) const = 0;
		virtual void PathConstraints(unsigned i, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const = 0;
		virtual void TerminalConstraints(const Eigen::VectorXd& x, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const = 0;
		virtual void Integrate(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x_next, Eigen::MatrixXd& A, Eigen::MatrixXd& B) const = 0;

	private:
		// Initialized _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		void UpdateStage(unsigned i);

		VectorMap w(unsigned i);
		VectorConstMap w(unsigned i) const;

		const double _sampleTime;

		const unsigned _Nu;
		const unsigned _Nx;
		const unsigned _Nz;
		const unsigned _Nt;
		const unsigned _Nd;
		const unsigned _NdT;
		
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

		// Lower state limit.
		Eigen::VectorXd _xMin;
		
		// Upper state limit.
		Eigen::VectorXd _xMax;

		// Lower terminal state limit.
		Eigen::VectorXd _terminalXMin;

		// Upper terminal state limit.
		Eigen::VectorXd _terminalXMax;
		
		// Lower input limit
		Eigen::VectorXd _uMin;
		
		// Upper input limit
		Eigen::VectorXd _uMax;		
	};

	template<class _Problem>
	ModelPredictiveController<_Problem>::ModelPredictiveController(unsigned state_dim, unsigned input_dim, unsigned n_path_constr, unsigned n_term_constr, double sample_time, unsigned Nt) :
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

	template<class _Problem>
	typename ModelPredictiveController<_Problem>::VectorMap ModelPredictiveController<_Problem>::w(unsigned i)
	{
		if(!(i < _Nt + 1))
			throw std::out_of_range("ModelPredictiveController<_Problem>::w(): index is out of range");

		return VectorMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
	}

	template<class _Problem>
	typename ModelPredictiveController<_Problem>::VectorConstMap ModelPredictiveController<_Problem>::w(unsigned i) const
	{
		if (!(i < _Nt + 1))
			throw std::out_of_range("ModelPredictiveController<_Problem>::w(): index is out of range");

		return VectorConstMap(_w.data() + i * _Nz, i < _Nt ? _Nz : _Nx);
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
		VectorXd u0(_Nu);
		u0.fill(0.);

		for (unsigned i = 0; i < _Nt; ++i)
			w(i) << x0, u0;

		w(_Nt) = x0;

		// Initialize QP
		UpdateQP();
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::Solve()
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

	template<class _Problem>
	void ModelPredictiveController<_Problem>::PrepareForNext()
	{
		/** prepare QP for next solution */
		//qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */
		//qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

		// Shift working point
		std::copy_n(_w.data() + _Nz, (_Nt - 1) * _Nz + _Nx, _w.data());

		// Calculate new matrices.
		UpdateQP();
	}

	template<class _Problem>
	Eigen::VectorXd ModelPredictiveController<_Problem>::getWorkingU(unsigned i) const
	{
		if (!(i < _Nt))
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
	void ModelPredictiveController<_Problem>::setXMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setXMin(): val has a wrong size");

		_xMin = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setXMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setXMax(): val has a wrong size");

		_xMax = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setTerminalXMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setTerminalXMin(): val has a wrong size");

		_terminalXMin = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setTerminalXMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nX())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setTerminalXMax(): val has a wrong size");

		_terminalXMax = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setUMin(const Eigen::VectorXd& val)
	{
		if (val.size() != nU())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setUMin(): val has a wrong size");

		_uMin = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setUMax(const Eigen::VectorXd& val)
	{
		if (val.size() != nU())
			throw std::invalid_argument("ModelPredictiveController<_Problem>::setUMax(): val has a wrong size");

		_uMax = val;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::setQPCallback(const QPCallback& cb)
	{
		_QPCallback = cb;
	}

	template<class _Problem>
	void ModelPredictiveController<_Problem>::UpdateStage(unsigned i)
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
