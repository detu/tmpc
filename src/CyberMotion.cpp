#pragma once

#include "CyberMotion.hpp"

#include <array>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ColMajorMatrix;

namespace CMS_output_jacobian
{
	void evaluate(const double * x, const double * u, double * y, double * dy_dx, double * dy_du);
	int getSparsity(int i, int *nrow, int *ncol, int **colind, int **row);
}

typedef std::array<int, 2> matrix_size_t;
matrix_size_t getSize(int i)
{
	matrix_size_t sz;
	int * colind, *row;
	int code = CMS_output_jacobian::getSparsity(i, &sz[0], &sz[1], &colind, &row);
	assert(code == 0);

	return sz;
}

CyberMotion::CyberMotion()
	: MotionPlatform(8)
{
}

void CyberMotion::Output(const double * x, const double * u, double * y, double * C, double * D) const
{
	// Call CasADi-generated code.
	CMS_output_jacobian::evaluate(x, u, y, C, D);
}

void CyberMotion::getDefaultAxesPosition(double * q) const
{
	// The "agile" position.
	const std::array<double, 8> q0 = { 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991 };
	std::copy(q0.begin(), q0.end(), q);
}

void CyberMotion::getAxesLimits(double * q_min, double * q_max, double * v_min, double * v_max, double * u_min, double * u_max) const
{
	const double data[8][6] = {
		{ 0.2, 9.3, -1.47, 1.47, -1.0780, 1.0780 },
		{ -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), -1.1802, 1.1802, -1.6762, 1.6762 },
		{ -2.1447, -0.8727, -0.9749, 0.9749, -1.1973, 1.1973 },
		{ -0.7540, 1.5415, -1.1802, 1.1802, -2.1893, 2.1893 },
		{ -3.0159, 3.0159, -1.2999, 1.2999, -0.5644, 0.5644 },
		{ -0.8727, 0.8727, -1.2999, 1.2999, -1.6249, 1.6249 },
		{ -3.0159, 3.0159, -2.0525, 2.0525, -1.3170, 1.3170 },
		// Full-range position for cabin axis is { 0.1938, 1.4839, -0.2450, 0.2450, -0.9800, 0.9800 }.
		// Below, the cabin axis position is limited to both rollers on the curve.
		// This is the only case implemented in current CasADi-generated code.
		{ 0.5700, 1.1603, -0.2450, 0.2450, -0.9800, 0.9800 }   
	};

	for (int i = 0; i < 8; ++i)
	{
		q_min[i] = data[i][0];	q_max[i] = data[i][1];
		v_min[i] = data[i][2];	v_max[i] = data[i][3];
		u_min[i] = data[i][4];	u_max[i] = data[i][5];
	}	
}