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
		MPC_Controller(const std::shared_ptr<MotionPlatform>& platform, double sample_time, unsigned Nt, const qpOptions_t * qpOptions = nullptr);
		~MPC_Controller();

		void InitWorkingPoint();
		void PrintQP(std::ostream& os) const;
		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		unsigned getNumberOfIntervals() const { return _Nt; }

	private:
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
		typedef Eigen::Map<RowMajorMatrix> RowMajorMatrixMap;
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;

		void InitQP();

		Eigen::MatrixXd getStateSpaceA() const;
		Eigen::MatrixXd getStateSpaceB() const;

		RowMajorMatrixMap H(unsigned i);
		VectorMap g(unsigned i);
		RowMajorMatrixMap C(unsigned i);
		VectorMap c(unsigned i);
		VectorMap zMin(unsigned i);
		VectorMap zMax(unsigned i);
		VectorMap xMin(unsigned i);
		VectorMap xMax(unsigned i);
		
		VectorMap z(unsigned i);
		VectorMap x(unsigned i);
		VectorMap u(unsigned i);
		VectorMap yRef(unsigned i);
		
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

		// Working point.
		// _z stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _z;

		// Reference output.
		// _yRef stores _Nt vectors of size _Ny
		std::vector<double> _yRef;
	};
}
