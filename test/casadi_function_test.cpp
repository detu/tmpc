#include <casadi_interface/GeneratedFunction.hpp>

#include "casadi_function_test_generated.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Dense>

#include "gtest_tools_eigen.hpp"

#include <tuple>

class CasADiFunctionTest : public ::testing::Test
{
protected:
	casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(f), 3, 2> fun_;
};

TEST_F(CasADiFunctionTest, incorrect_n_inputs_throws)
{
	typedef casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(f), 4, 2> wrong_fun_type;
	ASSERT_THROW(wrong_fun_type(), std::logic_error);
}

TEST_F(CasADiFunctionTest, incorrect_n_outputs_throws)
{
	typedef casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(f), 3, 3> wrong_fun_type;
	ASSERT_THROW(wrong_fun_type(), std::logic_error);
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
template <typename Function, typename... Args>
void call(Function const& f, Args&&... args)
{
	auto args_ = std::forward_as_tuple(std::forward<Args>(args)...);

	std::index_sequence<0, 1, 2> ind_in;
	std::index_sequence<0, 1> ind_out;
}

TEST_F(CasADiFunctionTest, matrix_argument_call_correct)
{
	Eigen::Matrix<double, 3, 2> const A = Eigen::Matrix<double, 3, 2>::Zero();
	Eigen::Matrix<double, 2, 2> const B = Eigen::Matrix<double, 2, 2>::Zero();
	Eigen::Matrix<double, 1, 1> const x = Eigen::Matrix<double, 1, 1>::Zero();

	Eigen::Matrix<double, 3, 2> X;
	Eigen::Matrix<double, 1, 2> Y;

	call(fun_, A, B, x, X, Y);
}
