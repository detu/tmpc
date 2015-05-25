#pragma once

#include <MotionPlatform.hpp>

#include <qpDUNES.h>

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
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<const RowMajorMatrix> RowMajorMatrixConstMap;

		MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt, const qpOptions_t * qpOptions = nullptr);
		~MPC_Controller();

		void InitWorkingPoint();
		void Solve(const double * x0, const double * y_ref);
		void getWorkingU(unsigned i, double * pu) const;
		void UpdateWorkingPoint();

		void PrintQP(std::ostream& os) const;
		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }
		unsigned getNumberOfIntervals() const { return _Nt; }
		RowMajorMatrixMap W(unsigned i);

	private:
		// Initialized _G, _y, _C, _c, _zMin, _zMax based on current working point _w.
		// Does not initialize g.
		void UpdateQP();

		void Integrate(const double * x, const double * u, double * x_next, double * A, double * B) const;
		Eigen::MatrixXd getStateSpaceA() const;
		Eigen::MatrixXd getStateSpaceB() const;

		RowMajorMatrixMap G(unsigned i);
		RowMajorMatrixMap H(unsigned i);
		VectorMap g(unsigned i);
		RowMajorMatrixMap C(unsigned i);
		VectorMap c(unsigned i);
		VectorMap zMin(unsigned i);
		VectorMap zMax(unsigned i);
		VectorMap xMin(unsigned i);
		VectorMap xMax(unsigned i);
		
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

		qpData_t _qpData;
		double _levenbergMarquardt;
		
		// Output weighting matrix
		// _W stores _Nt matrices of size _Ny x _Ny
		std::vector<double> _W;

		// _G stores _Nt row_major matrices of size _Ny x _Nz
		std::vector<double> _G;
		
		// _H stores _Nt row-major matrices of size _Nz x _Nz and 1 matrix of size _Nx x _Nx.
		std::vector<double> _H;

		// _g stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _g;

		// _C stores _Nt row-major matrices of size _Nx x _Nz
		std::vector<double> _C;

		// _c stores _Nt vectors of size _Nx
		std::vector<double> _c;

		// _zMin stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMin;

		// _zMax stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zMax;

		// Primal optimal solution.
		// _zOpt stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _zOpt;

		// Working point (linearization point).
		// _w stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _w;

		// Output at working point.
		// _y stores _Nt vectors of size _Ny.
		// Important: _y is column-major.
		Eigen::MatrixXd _y;

		// Reference output.
		// _yRef stores _Nt vectors of size _Ny
		//std::vector<double> _yRef;
	};
}
