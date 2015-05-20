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

class MotionLimits
{
public:
	MotionLimits(double q_min, double q_max, double v_min, double v_max, double u_min, double u_max)
		: _qMin(q_min), _qMax(q_max), _vMin(v_min), _vMax(v_max), _uMin(u_min), _uMax(u_max)
	{
		if (!(q_min < q_max && v_min < v_max && u_min < u_max))
			throw std::invalid_argument("MotionLimits::MotionLimits(): invalid motion limits (min must be less than max)");
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
	MotionPlatform(const std::vector<MotionLimits> limits)
		: System(static_cast<unsigned>(limits.size()), static_cast<unsigned>(2 * limits.size()), 6)
		, _gravity(0., 0., -9.81)
	{
	}

	unsigned getNumberOfAxes() const
	{
		return static_cast<unsigned>(_axesLimits.size());
	}

	const std::vector<MotionLimits>& getAxesLimits() const
	{
		return _axesLimits;
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

	virtual void Output(const ConstVectorMap& x, const ConstVectorMap& u, VectorMap& y, MatrixMap& C, MatrixMap& D) const = 0;

private:
	std::vector<MotionLimits> _axesLimits;
	Eigen::Vector3d _gravity;
};