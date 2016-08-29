#include <casadi_interface/GeneratedFunction.hpp>

#include "casadi_function_test_generated.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Dense>

#include "gtest_tools_eigen.hpp"

class CasADiFunctionTest : public ::testing::Test
{
protected:
	casadi_interface::GeneratedFunction fun_ {CASADI_GENERATED_FUNCTION_INTERFACE(f)};
};

TEST_F(CasADiFunctionTest, n_in_correct)
{
	EXPECT_EQ(fun_.n_in(), 3);
}

TEST_F(CasADiFunctionTest, n_out_correct)
{
	EXPECT_EQ(fun_.n_out(), 2);
}

TEST_F(CasADiFunctionTest, n_row_in_correct)
{
	EXPECT_EQ(fun_.n_row_in(0), 3);
	EXPECT_EQ(fun_.n_row_in(1), 2);
	EXPECT_EQ(fun_.n_row_in(2), 1);
}

TEST_F(CasADiFunctionTest, n_col_in_correct)
{
	EXPECT_EQ(fun_.n_col_in(0), 2);
	EXPECT_EQ(fun_.n_col_in(1), 2);
	EXPECT_EQ(fun_.n_col_in(2), 1);
}

TEST_F(CasADiFunctionTest, n_row_out_correct)
{
	EXPECT_EQ(fun_.n_row_out(0), 3);
	EXPECT_EQ(fun_.n_row_out(1), 1);
}

TEST_F(CasADiFunctionTest, n_col_out_correct)
{
	EXPECT_EQ(fun_.n_col_out(0), 2);
	EXPECT_EQ(fun_.n_col_out(1), 2);
}

TEST_F(CasADiFunctionTest, pointer_argument_call_correct)
{
	Eigen::Matrix<double, 3, 2> const A = Eigen::Matrix<double, 3, 2>::Zero();
	Eigen::Matrix<double, 2, 2> const B = Eigen::Matrix<double, 2, 2>::Zero();
	double const x = 0.;

	Eigen::Matrix<double, 3, 2> X;
	Eigen::Matrix<double, 1, 2> Y;

	fun_({A.data(), B.data(), &x}, {X.data(), Y.data()});
}

/*
TEST_F(CasADiFunctionTest, matrix_argument_call_correct)
{
	Eigen::Matrix<double, 3, 2> const A = Eigen::Matrix<double, 3, 2>::Zero();
	Eigen::Matrix<double, 2, 2> const B = Eigen::Matrix<double, 2, 2>::Zero();
	Eigen::Matrix<double, 1, 1> const x = Eigen::Matrix<double, 1, 1>::Zero();

	Eigen::Matrix<double, 3, 2> X;
	Eigen::Matrix<double, 1, 2> Y;

	fun_(A, B, x, X, Y);
}
*/
