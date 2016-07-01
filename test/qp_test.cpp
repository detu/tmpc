/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "../include/qp/MultiStageQuadraticProblem.hpp"
#include "../include/qp/HPMPCProblem.hpp"
#include "../include/qp/Printing.hpp"

#include <gtest/gtest.h>

#include <iostream>

template <typename Matrix>
bool isZero(Eigen::MatrixBase<Matrix> const& m)
{
	return m == Matrix::Zero();
}

template <typename Matrix>
Matrix random()
{
	Matrix m;

	for (std::size_t i = 0; i < m.rows(); ++i)
		for (std::size_t j = 0; j < m.cols(); ++j)
			m(i, j) = rand();

	return m;
}


template <typename QP>
class MultiStageQuadraticProblemTest : public ::testing::Test
{
public:
	typedef QP Problem;
	//typedef typename Problem::StateInputVector StateInputVector;
	//typedef typename Problem::StageHessianMatrix StageHessianMatrix;
	//typedef typename Problem::InterStageMatrix InterStageMatrix;

protected:
	/*
	unsigned const NX = n_x<QP>();
	unsigned const NU = n_u<QP>();
	unsigned const NC = n_d<QP>();
	unsigned const NCT = n_d_end<QP>();
	*/
	typedef Eigen::Matrix<double, QP::NX, 1> StateVector;
	typedef Eigen::Matrix<double, QP::NU, 1> InputVector;
	typedef Eigen::Matrix<double, QP::NX, QP::NX> StateStateMatrix;
	typedef Eigen::Matrix<double, QP::NX, QP::NU> StateInputMatrix;
	typedef Eigen::Matrix<double, QP::NU, QP::NU> InputInputMatrix;

	unsigned const NT = 2;

	MultiStageQuadraticProblemTest()
	:	qp(NT)
	{
	}

	Problem qp;
};

typedef ::testing::Types<
		tmpc::MultiStageQuadraticProblem<2, 1, 0, 0>
,		tmpc::HPMPCProblem              <2, 1, 0, 0>
	> QPTypes;

TYPED_TEST_CASE(MultiStageQuadraticProblemTest, QPTypes);

TYPED_TEST(MultiStageQuadraticProblemTest, get_set_interface_works)
{
	std::vector<typename TestFixture::StateStateMatrix> Q(this->NT + 1);
	std::vector<typename TestFixture::StateVector> x_min(this->NT + 1), x_max(this->NT + 1);

	std::vector<typename TestFixture::StateInputMatrix> S(this->NT + 1);
	std::vector<typename TestFixture::InputInputMatrix> R(this->NT + 1);
	std::vector<typename TestFixture::StateStateMatrix> A(this->NT);
	std::vector<typename TestFixture::StateInputMatrix> B(this->NT);
	std::vector<typename TestFixture::StateVector> b(this->NT);
	std::vector<typename TestFixture::InputVector> u_min(this->NT), u_max(this->NT);

	// Writing random data
	for (std::size_t i = 0; i <= this->NT; ++i)
	{
		this->qp.set_Q(i, Q[i] = random<typename TestFixture::StateStateMatrix>());
		this->qp.set_x_min(i, x_min[i] = random<typename TestFixture::StateVector>());
		this->qp.set_x_max(i, x_max[i] = random<typename TestFixture::StateVector>());
	}

	for (std::size_t i = 0; i < this->NT; ++i)
	{
		this->qp.set_S(i, S[i] = random<typename TestFixture::StateInputMatrix>());
		this->qp.set_R(i, R[i] = random<typename TestFixture::InputInputMatrix>());
		this->qp.set_A(i, A[i] = random<typename TestFixture::StateStateMatrix>());
		this->qp.set_B(i, B[i] = random<typename TestFixture::StateInputMatrix>());
		this->qp.set_b(i, b[i] = random<typename TestFixture::StateVector>());
		this->qp.set_u_min(i, u_min[i] = random<typename TestFixture::InputVector>());
		this->qp.set_u_max(i, u_max[i] = random<typename TestFixture::InputVector>());
	}

	// Reading the data and checking that they are the same that we wrote
	for (std::size_t i = 0; i <= this->NT; ++i)
	{
		EXPECT_EQ(this->qp.get_Q(i), Q[i]);
		EXPECT_EQ(this->qp.get_x_min(i), x_min[i]);
		EXPECT_EQ(this->qp.get_x_max(i), x_max[i]);
	}

	for (std::size_t i = 0; i < this->NT; ++i)
	{
		EXPECT_EQ(this->qp.get_S(i), S[i]);
		EXPECT_EQ(this->qp.get_R(i), R[i]);
		EXPECT_EQ(this->qp.get_A(i), A[i]);
		EXPECT_EQ(this->qp.get_B(i), B[i]);
		EXPECT_EQ(this->qp.get_b(i), b[i]);
		EXPECT_EQ(this->qp.get_u_min(i), u_min[i]);
		EXPECT_EQ(this->qp.get_u_max(i), u_max[i]);
	}
}

TYPED_TEST(MultiStageQuadraticProblemTest, qp_interface_works)
{
	auto const NZ = TestFixture::Problem::NX + TestFixture::Problem::NU;

	typedef Eigen::Matrix<double, NZ, 1> StateInputVector;
	typedef Eigen::Matrix<double, TestFixture::Problem::NX, 1> StateVector;
	typedef Eigen::Matrix<double, TestFixture::Problem::NU, 1> InputVector;
	typedef Eigen::Matrix<double, NZ, NZ> StageHessianMatrix;
	typedef Eigen::Matrix<double, TestFixture::Problem::NX, NZ> InterStageMatrix;

	set_xu_min(this->qp, 0, -1.);	set_xu_max(this->qp, 0, 1.);
	set_xu_min(this->qp, 1, -1.);	set_xu_max(this->qp, 1, 1.);
	set_x_end_min(this->qp, -1.);	set_x_end_max(this->qp, 1.);

	//std::cout << "******** QP *********" << std::endl;
	//Print_MATLAB(std::cout, qp, "qp");

	EXPECT_EQ(get_xu_min(this->qp, 0), StateInputVector::Constant(-1.));
	EXPECT_EQ(get_xu_max(this->qp, 0), StateInputVector::Constant( 1.));
	EXPECT_EQ(get_xu_min(this->qp, 1), StateInputVector::Constant(-1.));
	EXPECT_EQ(get_xu_max(this->qp, 1), StateInputVector::Constant( 1.));
	EXPECT_EQ(get_x_end_min(this->qp), StateVector::Constant(-1.));
	EXPECT_EQ(get_x_end_max(this->qp), StateVector::Constant( 1.));

	// Stage 0
	StageHessianMatrix H0;
	H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H0 = H0.transpose() * H0;	// Make positive definite.

	const Eigen::MatrixXd Q0 = H0.topLeftCorner(this->qp.nX(), this->qp.nX());
	const Eigen::MatrixXd R0 = H0.bottomRightCorner(this->qp.nU(), this->qp.nU());
	const Eigen::MatrixXd S0 = H0.topRightCorner(this->qp.nX(), this->qp.nU());
	const Eigen::MatrixXd S0T = H0.bottomLeftCorner(this->qp.nU(), this->qp.nX());

	Eigen::MatrixXd A0(this->qp.nX(), this->qp.nX());
	A0 << 1, 1, 0, 1;

	Eigen::MatrixXd B0(this->qp.nX(), this->qp.nU());
	B0 << 0.5, 1.0;

	Eigen::VectorXd a0(this->qp.nX());
	a0 << 1, 2;

	// Stage 1
	StageHessianMatrix H1;
	H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H1 = H1.transpose() * H1;	// Make positive definite.

	const Eigen::MatrixXd Q1 = H1.topLeftCorner(this->qp.nX(), this->qp.nX());
	const Eigen::MatrixXd R1 = H1.bottomRightCorner(this->qp.nU(), this->qp.nU());
	const Eigen::MatrixXd S1 = H1.topRightCorner(this->qp.nX(), this->qp.nU());
	const Eigen::MatrixXd S1T = H1.bottomLeftCorner(this->qp.nU(), this->qp.nX());

	Eigen::MatrixXd A1(this->qp.nX(), this->qp.nX());
	A1 << 1, 1, 0, 1;

	Eigen::MatrixXd B1(this->qp.nX(), this->qp.nU());
	B1 << 0.5, 1.0;

	Eigen::VectorXd a1(this->qp.nX());
	a1 << 1, 2;

	// Stage 2
	Eigen::MatrixXd H2(this->qp.nX(), this->qp.nX());
	H2 << 1, 2, 3, 4;
	H2 = H2.transpose() * H2;	// Make positive definite.

	const Eigen::MatrixXd Q2 = H2.topLeftCorner(this->qp.nX(), this->qp.nX());

	// Setup QP
	set_H(this->qp, 0, H0);
	EXPECT_EQ(get_H(this->qp, 0), H0);

	set_H(this->qp, 1, H1);
	EXPECT_EQ(get_H(this->qp, 1), H1);

	set_Q_end(this->qp, H2);
	EXPECT_EQ(get_Q_end(this->qp), H2);

	InterStageMatrix C0;
	C0 << A0, B0;
	set_AB(this->qp, 0, C0);
	EXPECT_EQ(get_AB(this->qp, std::size_t(0)), C0);
	this->qp.set_b(0, a0);
	EXPECT_EQ(this->qp.get_b(0), a0);

	InterStageMatrix C1;
	C1 << A1, B1;
	set_AB(this->qp, 1, C1);
	EXPECT_EQ(get_AB(this->qp, 1), C1);

	this->qp.set_b(1, a1);
	EXPECT_EQ(this->qp.get_b(1), a1);
	//EXPECT_EQ(Eigen::Map<Problem::StateVector const>(this->qp.b_data()[1]), a1);
}
