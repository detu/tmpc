#pragma once

#include <MultiStageQP.hpp>
#include <CondensingSolver.hpp>
#include <MotionPlatform.hpp>

#include <Eigen/Dense>

#include <vector>
#include <memory>
#include <ostream>

namespace rtmc
{
	class MPC_Controller
	{
	public:
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt);
		~MPC_Controller();

		void InitWorkingPoint();
		void Solve();

		void EmbedInitialValue(const double * px0);

		void SetReference(const double * py_ref);

		void getWorkingU(unsigned i, double * pu) const;
		void PrepareForNext();

		void PrintQP_C(std::ostream& os) const;
		void PrintQP_zMax_C(std::ostream& log_stream) const;
		void PrintQP_zMin_C(std::ostream& log_stream) const;
		void PrintQP_MATLAB(std::ostream& log_stream) const;
		
		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		double getWashoutFactor() const { return _washoutFactor; }
		void setWashoutFactor(double val) { _washoutFactor = val; }

		unsigned getNumberOfIntervals() const { return _Nt; }
		double getSampleTime() const;

		RowMajorMatrixMap W(unsigned i);

	private:
		// Initialized _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		Eigen::VectorXd getWashoutState() const;

		void Integrate(const double * x, const double * u, double * x_next, double * A, double * B) const;
		Eigen::MatrixXd getStateSpaceA() const;
		Eigen::MatrixXd getStateSpaceB() const;

		RowMajorMatrixMap G(unsigned i);		
		VectorMap w(unsigned i);
		//VectorMap yRef(unsigned i);
		
		const std::shared_ptr<MotionPlatform> _platform;
		double _sampleTime;

		unsigned _Nq;
		unsigned _Nu;
		unsigned _Nx;
		unsigned _Ny;
		unsigned _Nz;
		unsigned _Nt;
		
		camels::MultiStageQP _QP;
		camels::CondensingSolver _Solver;

		double _levenbergMarquardt;

		// The more the washout factor, the more penalty for the terminal state to be far from the default (washout) position.
		double _washoutFactor = 0.;

		// Output weighting matrix
		// _W stores _Nt matrices of size _Ny x _Ny
		std::vector<double> _W;

		// _G stores _Nt row_major matrices of size _Ny x _Nz
		std::vector<double> _G;
		
		// Primal optimal solution.
		// _zOpt stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zOpt;

		// Working point (linearization point).
		// _w stores _Nt vectors of size _Nz and 1 vector of size _Nx
		Eigen::VectorXd _w;

		// Output at working point.
		// _y stores _Nt vectors of size _Ny.
		// Important: _y is column-major.
		Eigen::MatrixXd _y;

		// Reference output.
		// _yRef stores _Nt vectors of size _Ny
		//std::vector<double> _yRef;
	};
}
