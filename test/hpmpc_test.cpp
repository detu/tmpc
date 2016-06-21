/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "../include/qp/HPMPCSolver.hpp"

#include <gtest/gtest.h>

#include <iostream>

unsigned const NX = 2;
unsigned const NU = 1;
unsigned const NC = 0;
unsigned const NCT = 0;
unsigned const NT = 2;

typedef tmpc::HPMPCSolver<NX, NU, NC, NCT> Solver;
typedef Solver::Problem Problem;
typedef Solver::Solution Solution;

namespace
{
	std::ostream& operator<<(std::ostream& os, Solution const& point)
	{
		//typedef typename camels::CondensingSolver<NX_, NU_, NC_, NCT_>::size_type size_type;
		typedef unsigned size_type;
		for (size_type i = 0; i < point.nT(); ++i)
			os << point.get_x(i) << "\t" << point.get_u(i) << std::endl;

		return os << point.get_xend() << std::endl;
	}
}

TEST(hpmpc_test, problem_test)
{
	Problem qp(NT);
	set_zMin(qp, 0, -1.);	set_zMax(qp, 0, 1.);
	set_zMin(qp, 1, -1.);	set_zMax(qp, 1, 1.);
	set_zendMin(qp, -1.);	set_zendMax(qp, 1.);

	EXPECT_EQ(get_zMin(qp, 0), Problem::StateInputVector::Constant(-1.));
	EXPECT_EQ(get_zMax(qp, 0), Problem::StateInputVector::Constant( 1.));
	EXPECT_EQ(get_zMin(qp, 1), Problem::StateInputVector::Constant(-1.));
	EXPECT_EQ(get_zMax(qp, 1), Problem::StateInputVector::Constant( 1.));
	EXPECT_EQ(get_zendMin(qp), Problem::StateVector::Constant(-1.));
	EXPECT_EQ(get_zendMax(qp), Problem::StateVector::Constant( 1.));

	EXPECT_EQ(Eigen::Map<Problem::StateInputVector const>(qp.lb_data()[0]), Problem::StateInputVector::Constant(-1.));
	EXPECT_EQ(Eigen::Map<Problem::StateInputVector const>(qp.ub_data()[0]), Problem::StateInputVector::Constant( 1.));
	EXPECT_EQ(Eigen::Map<Problem::StateInputVector const>(qp.lb_data()[1]), Problem::StateInputVector::Constant(-1.));
	EXPECT_EQ(Eigen::Map<Problem::StateInputVector const>(qp.ub_data()[1]), Problem::StateInputVector::Constant( 1.));
	EXPECT_EQ(Eigen::Map<Problem::StateVector const>(qp.lb_data()[2]), Problem::StateVector::Constant(-1.));
	EXPECT_EQ(Eigen::Map<Problem::StateVector const>(qp.ub_data()[2]), Problem::StateVector::Constant( 1.));

	// Stage 0
	Problem::StageHessianMatrix H0;
	H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H0 = H0.transpose() * H0;	// Make positive definite.

	const Eigen::MatrixXd Q0 = H0.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R0 = H0.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S0 = H0.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S0T = H0.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A0(qp.nX(), qp.nX());
	A0 << 1, 1, 0, 1;

	Eigen::MatrixXd B0(qp.nX(), qp.nU());
	B0 << 0.5, 1.0;

	Eigen::VectorXd a0(qp.nX());
	a0 << 1, 2;

	// Stage 1
	Problem::StageHessianMatrix H1;
	H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H1 = H1.transpose() * H1;	// Make positive definite.

	const Eigen::MatrixXd Q1 = H1.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R1 = H1.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S1 = H1.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S1T = H1.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A1(qp.nX(), qp.nX());
	A1 << 1, 1, 0, 1;

	Eigen::MatrixXd B1(qp.nX(), qp.nU());
	B1 << 0.5, 1.0;

	Eigen::VectorXd a1(qp.nX());
	a1 << 1, 2;

	// Stage 2
	Eigen::MatrixXd H2(qp.nX(), qp.nX());
	H2 << 1, 2, 3, 4;
	H2 = H2.transpose() * H2;	// Make positive definite.

	const Eigen::MatrixXd Q2 = H2.topLeftCorner(qp.nX(), qp.nX());

	// Setup QP
	set_H(qp, 0, H0);
	EXPECT_EQ(get_H(qp, 0), H0);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_QMatrix const>(qp.Q_data()[0]), Q0);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_RMatrix const>(qp.R_data()[0]), R0);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_SMatrix const>(qp.S_data()[0]), S0.transpose());

	set_H(qp, 1, H1);
	EXPECT_EQ(get_H(qp, 1), H1);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_QMatrix const>(qp.Q_data()[1]), Q1);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_RMatrix const>(qp.R_data()[1]), R1);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_SMatrix const>(qp.S_data()[1]), S1.transpose());

	set_Hend(qp, H2);
	EXPECT_EQ(get_Hend(qp), H2);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_QMatrix const>(qp.Q_data()[2]), Q2);
	EXPECT_EQ(qp.R_data()[2], nullptr);
	EXPECT_EQ(qp.S_data()[2], nullptr);

	EXPECT_NE(qp.q_data()[2], nullptr);
	EXPECT_EQ(qp.r_data()[2], nullptr);

	Problem::InterStageMatrix C0;
	C0 << A0, B0;
	set_C(qp, 0, C0);
	EXPECT_EQ(get_C(qp, std::size_t(0)), C0);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_AMatrix const>(qp.A_data()[0]), A0);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_BMatrix const>(qp.B_data()[0]), B0);

	set_c(qp, 0, a0);
	EXPECT_EQ(get_c(qp, 0), a0);
	EXPECT_EQ(Eigen::Map<Problem::StateVector const>(qp.b_data()[0]), a0);

	Problem::InterStageMatrix C1;
	C1 << A1, B1;
	set_C(qp, 1, C1);
	EXPECT_EQ(get_C(qp, 1), C1);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_AMatrix const>(qp.A_data()[1]), A1);
	EXPECT_EQ(Eigen::Map<Problem::HPMPC_BMatrix const>(qp.B_data()[1]), B1);

	set_c(qp, 1, a1);
	EXPECT_EQ(get_c(qp, 1), a1);
	EXPECT_EQ(Eigen::Map<Problem::StateVector const>(qp.b_data()[1]), a1);
 }

TEST(hpmpc_test, solve_test)
{
	Problem qp(NT);
	set_zMin(qp, 0, -1.);	set_zMax(qp, 0, 1.);
	set_zMin(qp, 1, -1.);	set_zMax(qp, 1, 1.);
	set_zendMin(qp, -1.);	set_zendMax(qp, 1.);

	// Stage 0
	Problem::StageHessianMatrix H0;
	H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H0 = H0.transpose() * H0;	// Make positive definite.

	const Eigen::MatrixXd Q0 = H0.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R0 = H0.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S0 = H0.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S0T = H0.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A0(qp.nX(), qp.nX());
	A0 << 1, 1, 0, 1;

	Eigen::MatrixXd B0(qp.nX(), qp.nU());
	B0 << 0.5, 1.0;

	Eigen::VectorXd a0(qp.nX());
	a0 << 1, 2;

	// Stage 1
	Problem::StageHessianMatrix H1;
	H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H1 = H1.transpose() * H1;	// Make positive definite.

	const Eigen::MatrixXd Q1 = H1.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R1 = H1.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S1 = H1.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S1T = H1.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A1(qp.nX(), qp.nX());
	A1 << 1, 1, 0, 1;

	Eigen::MatrixXd B1(qp.nX(), qp.nU());
	B1 << 0.5, 1.0;

	Eigen::VectorXd a1(qp.nX());
	a1 << 1, 2;

	// Stage 2
	Eigen::MatrixXd H2(qp.nX(), qp.nX());
	H2 << 1, 2, 3, 4;
	H2 = H2.transpose() * H2;	// Make positive definite.

	const Eigen::MatrixXd Q2 = H2.topLeftCorner(qp.nX(), qp.nX());

	// Setup QP
	set_H(qp, 0, H0);
	set_H(qp, 1, H1);
	set_Hend(qp, H2);

	Problem::InterStageMatrix C0;
	C0 << A0, B0;
	set_C(qp, 0, C0);
	set_c(qp, 0, a0);

	Problem::InterStageMatrix C1;
	C1 << A1, B1;
	set_C(qp, 1, C1);
	set_c(qp, 1, a1);

	Solver solver(qp.nT());
	Solution solution(NT);
	solver.Solve(qp, solution);

	std::cout << "-- sol (multistage) --" << std::endl << solution << std::endl;
 }
