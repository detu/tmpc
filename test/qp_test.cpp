/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "gtest_tools_eigen.hpp"

#include <tmpc/qp/HPMPCProblem.hpp>
#include <tmpc/qp/Printing.hpp>
#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/QuadraticProblem.hpp>

#include <gtest/gtest.h>

#include <iostream>
#include <array>

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

template <unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
class FixedSizeQpOasesProblem : public tmpc::QpOasesProblem
{
public:
	static auto constexpr NX = NX_;
	static auto constexpr NU = NU_;
	static auto constexpr NC = NC_;
	static auto constexpr NCT = NCT_;

	/*
	static unsigned constexpr NX = NX_;
	static unsigned const NU = NU_;
	static unsigned const NC = NC_;
	static unsigned const NCT = NCT_;
	*/

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;

	FixedSizeQpOasesProblem(std::size_t N)
	:	tmpc::QpOasesProblem(Sizes(N))
	{
	}

	static unsigned constexpr nX() { return NX; }
	static unsigned constexpr nU() { return NU; }

private:
	std::vector<tmpc::QpSize> Sizes(std::size_t N)
	{
		std::vector<tmpc::QpSize> sz;
		sz.reserve(N + 1);

		std::fill_n(std::back_inserter(sz), N, tmpc::QpSize(NX, NU, 0));
		sz.emplace_back(tmpc::QpSize::size_type{NX}, 0, 0);
		//sz.push_back(tmpc::QpSize(NX, 0, 0));

		return sz;
	}
};

template <typename QP>
class QuadraticProblemTest : public ::testing::Test
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

	QuadraticProblemTest()
	:	qp(NT)
	{
	}

	Problem qp;
};

typedef ::testing::Types<
	tmpc::QuadraticProblemEigen<2, 1, 0, 0>,
	tmpc::HPMPCProblem<2, 1, 0, 0>,
	FixedSizeQpOasesProblem<2, 1, 0, 0>
	> QPTypes;

TYPED_TEST_CASE(QuadraticProblemTest, QPTypes);

TYPED_TEST(QuadraticProblemTest, get_set_interface_works)
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
		EXPECT_EQ(print_wrap(this->qp.get_Q(i)), print_wrap(Q[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(this->qp.get_x_min(i)), print_wrap(x_min[i]));
		EXPECT_EQ(print_wrap(this->qp.get_x_max(i)), print_wrap(x_max[i]));
	}

	for (std::size_t i = 0; i < this->NT; ++i)
	{
		EXPECT_EQ(print_wrap(this->qp.get_S(i)), print_wrap(S[i]));
		EXPECT_EQ(print_wrap(this->qp.get_R(i)), print_wrap(R[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(this->qp.get_A(i)), print_wrap(A[i]));
		EXPECT_EQ(print_wrap(this->qp.get_B(i)), print_wrap(B[i]));
		EXPECT_EQ(print_wrap(this->qp.get_b(i)), print_wrap(b[i]));
		EXPECT_EQ(print_wrap(this->qp.get_u_min(i)), print_wrap(u_min[i]));
		EXPECT_EQ(print_wrap(this->qp.get_u_max(i)), print_wrap(u_max[i]));
	}
}

TYPED_TEST(QuadraticProblemTest, qp_interface_works)
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

template <typename QP>
class QuadraticProblemTestBlaze : public ::testing::Test
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
	typedef typename Problem::Scalar Scalar;

	typedef blaze::DynamicVector<Scalar> Vector;
	typedef blaze::DynamicMatrix<Scalar> Matrix;

	QuadraticProblemTestBlaze()
	:	size_{
			tmpc::QpSize(2, 1, 0),
			tmpc::QpSize(2, 1, 0),
			tmpc::QpSize(2, 0, 0)
		}
	,	qp_(size_.begin(), size_.end())
	{
	}

	std::array<tmpc::QpSize, 3> const size_;
	Problem qp_;
};

typedef ::testing::Types<
	tmpc::QuadraticProblemBlaze<double>
	> QPTypesBlaze;

TYPED_TEST_CASE(QuadraticProblemTestBlaze, QPTypesBlaze);

TYPED_TEST(QuadraticProblemTestBlaze, get_set_interface_works)
{
	auto const N = this->size_.size();

	std::vector<typename TestFixture::Matrix> Q(N);
	std::vector<typename TestFixture::Vector> q(N);
	std::vector<typename TestFixture::Vector> r(N);

	std::vector<typename TestFixture::Matrix> S(N);
	std::vector<typename TestFixture::Matrix> R(N);
	std::vector<typename TestFixture::Matrix> A(N);
	std::vector<typename TestFixture::Matrix> B(N);
	std::vector<typename TestFixture::Vector> b(N);

	std::vector<typename TestFixture::Vector> x_min(N), x_max(N);
	std::vector<typename TestFixture::Vector> u_min(N), u_max(N);

	// Writing random data
	blaze::Rand<typename TestFixture::Matrix> rand_matrix;
	blaze::Rand<typename TestFixture::Vector> rand_vector;

	for (std::size_t i = 0; i < N; ++i)
	{
		auto const& sz = this->size_[i];
		auto const nx1 = i + 1 < N ? this->size_[i + 1].nx() : 0;
		auto& stage = this->qp_[i];

		stage.set_Q(Q[i] = rand_matrix.generate(sz.nx(), sz.nx()));
		stage.set_R(R[i] = rand_matrix.generate(sz.nu(), sz.nu()));
		stage.set_S(S[i] = rand_matrix.generate(sz.nx(), sz.nu()));
		stage.set_q(q[i] = rand_vector.generate(sz.nx()));
		stage.set_r(r[i] = rand_vector.generate(sz.nx()));

		stage.set_A(A[i] = rand_matrix.generate(nx1, sz.nx()));
		stage.set_B(B[i] = rand_matrix.generate(nx1, sz.nu()));
		stage.set_b(b[i] = rand_vector.generate(nx1));

		stage.set_lbx(x_min[i] = rand_vector.generate(sz.nx()));
		stage.set_ubx(x_max[i] = rand_vector.generate(sz.nx()));
		stage.set_lbu(u_min[i] = rand_vector.generate(sz.nu()));
		stage.set_ubu(u_max[i] = rand_vector.generate(sz.nu()));	}

	// Reading the data and checking that they are the same that we wrote
	for (std::size_t i = 0; i < N; ++i)
	{
		auto const& stage = this->qp_[i];

		EXPECT_EQ(print_wrap(stage.get_Q()), print_wrap(Q[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(stage.get_R()), print_wrap(R[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(stage.get_S()), print_wrap(S[i]));
		EXPECT_EQ(print_wrap(stage.get_q()), print_wrap(q[i]));
		EXPECT_EQ(print_wrap(stage.get_r()), print_wrap(r[i]));

		EXPECT_EQ(print_wrap(stage.get_A()), print_wrap(A[i]));
		EXPECT_EQ(print_wrap(stage.get_B()), print_wrap(B[i]));
		EXPECT_EQ(print_wrap(stage.get_b()), print_wrap(b[i]));

		EXPECT_EQ(print_wrap(stage.get_lbx()), print_wrap(x_min[i]));
		EXPECT_EQ(print_wrap(stage.get_ubx()), print_wrap(x_max[i]));
		EXPECT_EQ(print_wrap(stage.get_lbu()), print_wrap(u_min[i]));
		EXPECT_EQ(print_wrap(stage.get_ubu()), print_wrap(u_max[i]));
	}
}
