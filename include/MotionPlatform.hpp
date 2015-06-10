#pragma once

#include <Eigen/Dense>

#include <vector>

class System
{
public:
	unsigned getInputDim() const
	{
		return _inputDim;
	}

	unsigned getStateDim() const
	{
		return _stateDim;
	}

	unsigned getOutputDim() const
	{
		return _outputDim;
	}

protected:
	System(unsigned nu, unsigned nx, unsigned ny)
		: _inputDim(nu), _stateDim(nx), _outputDim(ny)
	{
	}

private:
	unsigned _inputDim;
	unsigned _stateDim;
	unsigned _outputDim;
};

namespace rtmc
{
	class MotionLimits
	{
	public:
		MotionLimits(double q_min, double q_max, double v_min, double v_max, double u_min, double u_max)
			: _qMin(q_min), _qMax(q_max), _vMin(v_min), _vMax(v_max), _uMin(u_min), _uMax(u_max)
		{
			if (!(q_min < q_max && v_min < v_max && u_min < u_max))
				throw std::invalid_argument("MotionLimits::MotionLimits(): invalid motion limits (min must be less than max)");
		}

		double getPositionMin() const
		{
			return _qMin;
		}

		double getPositionMax() const
		{
			return _qMax;
		}

		double getVelocityMin() const
		{
			return _vMin;
		}

		double getVelocityMax() const
		{
			return _vMax;
		}

		double getAccelerationMin() const
		{
			return _uMin;
		}

		double getAccelerationMax() const
		{
			return _uMax;
		}

	private:
		double _qMin;
		double _qMax;
		double _vMin;
		double _vMax;
		double _uMin;
		double _uMax;
	};

	class MotionPlatform : public System
	{
	public:
		MotionPlatform(unsigned Nq)
			: System(Nq, 2 * Nq, 6)
			, _gravity(0., 0., -9.81)
			, _numberOfAxes(Nq)
		{
		}

		unsigned getNumberOfAxes() const
		{
			return _numberOfAxes;
		}

		const Eigen::Vector3d& getGravity() const
		{
			return _gravity;
		}

		void setGravity(const Eigen::Vector3d& g)
		{
			_gravity = g;
		}

		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> DynamicStride;
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> ConstVectorMap;
		typedef Eigen::Map<Eigen::MatrixXd, 0, DynamicStride> MatrixMap;

		virtual void getAxesLimits(double * q_min, double * q_max, double * v_min, double * v_max, double * u_min, double * u_max) const = 0;
		virtual void getDefaultAxesPosition(double * q) const = 0;

		// Output matrices storage order is column-major.
		virtual void Output(const double * x, const double * u, double * y, double * C = nullptr, double * D = nullptr) const = 0;

	private:
		Eigen::Vector3d _gravity;
		const unsigned _numberOfAxes;
	};
}