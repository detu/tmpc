#pragma once

#include <MultiStageQP.hpp>
#include <CondensingSolver.hpp>

#include <Eigen/Dense>

#include <vector>
#include <memory>
#include <ostream>

namespace camels
{
	class MPC_Controller
	{
	public:
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		MPC_Controller(unsigned state_dim, unsigned input_dim, double sample_time, unsigned Nt);
		~MPC_Controller();

		void InitWorkingPoint(const Eigen::VectorXd& x0);
		void Solve();

		void EmbedInitialValue(const double * px0);

		void getWorkingU(unsigned i, double * pu) const;
		void PrepareForNext();

		void PrintQP_C(std::ostream& os) const;
		void PrintQP_zMax_C(std::ostream& log_stream) const;
		void PrintQP_zMin_C(std::ostream& log_stream) const;
		void PrintQP_MATLAB(std::ostream& log_stream) const;
		
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

	protected:
		virtual void LagrangeTerm(const Eigen::MatrixXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) = 0;
		virtual void MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) = 0;
		virtual void Integrate(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& x_next, Eigen::MatrixXd& A, Eigen::MatrixXd& B) const = 0;

	private:
		// Initialized _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		VectorMap w(unsigned i);

		double _sampleTime;

		unsigned _Nu;
		unsigned _Nx;
		unsigned _Nz;
		unsigned _Nt;
		
		camels::MultiStageQP _QP;
		camels::CondensingSolver _Solver;

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
}
