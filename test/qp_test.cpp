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

using namespace tmpc;

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
	typedef typename Problem::Scalar Scalar;

	typedef DynamicVector<Scalar> Vector;
	typedef DynamicMatrix<Scalar> Matrix;

	QuadraticProblemTest()
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
	QuadraticProblem<double>,
	HPMPCProblem<double>,
	QpOasesProblem //<double>
	> QPTypes;

TYPED_TEST_CASE(QuadraticProblemTest, QPTypes);

TYPED_TEST(QuadraticProblemTest, get_set_interface_works)
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
	Rand<typename TestFixture::Matrix> rand_matrix;
	Rand<typename TestFixture::Vector> rand_vector;

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
