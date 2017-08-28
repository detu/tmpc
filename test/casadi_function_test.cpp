#include <tmpc/casadi_interface/GeneratedFunction.hpp>
#include <tmpc/test_tools.hpp>

#include "casadi_function_test_generated.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Dense>

#include <tuple>

class CasADiFunctionTest : public ::testing::Test
{
protected:
	casadi_interface::GeneratedFunction fun_ {f_functions()};
};

TEST_F(CasADiFunctionTest, incorrect_n_inputs_throws)
{
	double x = 0.;
	ASSERT_THROW(fun_({&x, &x, &x, &x}, {&x, &x}), std::invalid_argument);
}

TEST_F(CasADiFunctionTest, incorrect_n_outputs_throws)
{
	double x = 0.;
	ASSERT_THROW(fun_({&x, &x, &x}, {&x, &x, &x}), std::invalid_argument);
}

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

// Some interesting ideas here: https://habrahabr.ru/post/228031/
template <typename Function, typename TupleIn, typename TupleOut>
void call(Function const& f, TupleIn&& in, TupleOut&& out)
{
	throw std::logic_error("Not implemented");
}

TEST_F(CasADiFunctionTest, DISABLED_matrix_argument_call_correct)
{
	Eigen::Matrix<double, 3, 2> const A = Eigen::Matrix<double, 3, 2>::Zero();
	Eigen::Matrix<double, 2, 2> const B = Eigen::Matrix<double, 2, 2>::Zero();
	Eigen::Matrix<double, 1, 1> const x = Eigen::Matrix<double, 1, 1>::Zero();

	Eigen::Matrix<double, 3, 2> X;
	Eigen::Matrix<double, 1, 2> Y;

	call(fun_, std::forward_as_tuple(A, B, x), std::forward_as_tuple(X, Y));
}
