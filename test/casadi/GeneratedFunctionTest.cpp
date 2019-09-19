#include <tmpc/casadi/GeneratedFunction.hpp>
#include <tmpc/Testing.hpp>

#include <test_functions.h>

#include <blaze/Math.h>

#include <tuple>
#include <stdexcept>
#include <random>


namespace tmpc :: testing
{
	using Real = double;
	using blaze::StaticMatrix;
	using blaze::StaticVector;
	using blaze::columnMajor;
	using blaze::rowVector;
	using blaze::columnVector;
	using namespace tmpc :: casadi;	


	class GeneratedFunctionTest 
	: 	public Test
	{
	protected:
		casadi::GeneratedFunction fun_ {f_functions()};
		using MatrixA = StaticMatrix<Real, 3, 2, columnMajor>;
		using MatrixB = StaticMatrix<Real, 2, 2, columnMajor>;
		using MatrixX = StaticMatrix<Real, 3, 2>;
		using VectorY = StaticVector<Real, 2, rowVector>;
	};


	TEST_F(GeneratedFunctionTest, incorrect_n_inputs_throws)
	{
		Real x = 0.;
		ASSERT_THROW(fun_({&x, &x, &x, &x}, {&x, &x}), std::invalid_argument);
	}


	TEST_F(GeneratedFunctionTest, incorrect_n_outputs_throws)
	{
		Real x = 0.;
		ASSERT_THROW(fun_({&x, &x, &x}, {&x, &x, &x}), std::invalid_argument);
	}


	TEST_F(GeneratedFunctionTest, n_in_correct)
	{
		EXPECT_EQ(fun_.n_in(), 3);
	}


	TEST_F(GeneratedFunctionTest, n_out_correct)
	{
		EXPECT_EQ(fun_.n_out(), 2);
	}


	TEST_F(GeneratedFunctionTest, n_row_in_correct)
	{
		EXPECT_EQ(fun_.n_row_in(0), 3);
		EXPECT_EQ(fun_.n_row_in(1), 2);
		EXPECT_EQ(fun_.n_row_in(2), 1);
	}


	TEST_F(GeneratedFunctionTest, n_col_in_correct)
	{
		EXPECT_EQ(fun_.n_col_in(0), 2);
		EXPECT_EQ(fun_.n_col_in(1), 2);
		EXPECT_EQ(fun_.n_col_in(2), 1);
	}


	TEST_F(GeneratedFunctionTest, n_row_out_correct)
	{
		EXPECT_EQ(fun_.n_row_out(0), 3);
		EXPECT_EQ(fun_.n_row_out(1), 1);
	}


	TEST_F(GeneratedFunctionTest, n_col_out_correct)
	{
		EXPECT_EQ(fun_.n_col_out(0), 2);
		EXPECT_EQ(fun_.n_col_out(1), 2);
	}


	TEST_F(GeneratedFunctionTest, testPointerArgumentCall)
	{
		Real const A[3 * 2] = {
			1., 3., 5.,
			2., 4., 6.
		};
		
		Real const B[2 * 2] = {
			7., 9.,
			8., 10.
		};

		Real x = 0.1;

		Real X[3 * 2];
		Real Y[2];

		fun_({A, B, &x}, {X, Y});

		EXPECT_THAT(X, ElementsAreArray({
			DoubleEq(2.5), DoubleEq(5.7), DoubleEq(8.9),
			DoubleEq(2.8), DoubleEq(6.4), DoubleEq(10.0)
		}));
		EXPECT_THAT(Y, ElementsAreArray({171., 192.}));
	}


	TEST_F(GeneratedFunctionTest, testCopyCtor)
	{
		StaticMatrix<Real, 3, 2, columnMajor> A {
			{1., 2.},
			{3., 4.},
			{5., 6.}
		};
		StaticMatrix<Real, 2, 2, columnMajor> B {
			{7., 8.},
			{9., 10.}
		};
		Real x = 0.1;

		//randomize(A);
		//randomize(B);
		//blaze::randomize(x);

		StaticMatrix<Real, 3, 2, columnMajor> X, X1;
		StaticVector<Real, 2, rowVector> Y, Y1;

		casadi::GeneratedFunction const fun_copy = fun_;
		fun_(std::tie(A, B, x), std::tie(X, Y));
		fun_copy(std::tie(A, B, x), std::tie(X1, Y1));

		TMPC_EXPECT_APPROX_EQ(X, (A * x) * B, 1e-12, 0.);
		TMPC_EXPECT_EQ(Y, (StaticVector<Real, 3, rowVector> {1., 1., 1.} * (A * B)));
		TMPC_EXPECT_APPROX_EQ(X1, (A * x) * B, 1e-12, 0.);
		TMPC_EXPECT_EQ(Y1, (StaticVector<Real, 3, rowVector> {1., 1., 1.} * (A * B)));
	}


	TEST_F(GeneratedFunctionTest, testOmpParallelEvaluation)
	{
		size_t const N = 5000;

		std::vector<MatrixA> A(N);
		std::vector<MatrixB> B(N);
		std::vector<Real> x(N);

		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dis(-1.0, 1.0);	

		for (size_t i = 0; i < N; ++i)
		{
			randomize(A[i]);
			randomize(B[i]);
			x[i] = dis(gen);
		}

		std::vector<MatrixX> X(N), X_ref(N);
		std::vector<VectorY> Y(N), Y_ref(N);

		// Evaluate the function sequentially.
		for (size_t i = 0; i < N; ++i)
			fun_({A[i].data(), B[i].data(), &x[i]}, {X_ref[i].data(), Y_ref[i].data()});

		// Evaluate the function in parallel. 
		// Set high number of threads to maximize probability of collisions.
		#pragma omp parallel num_threads(100) firstprivate(fun_)
		{
			#pragma omp for
			for (size_t i = 0; i < N; ++i)
				fun_({A[i].data(), B[i].data(), &x[i]}, {X[i].data(), Y[i].data()});
		}

		// Check if the results are the same.
		for (size_t i = 0; i < N; ++i)
		{
			TMPC_EXPECT_EQ(X[i], X_ref[i]) << "at i=" << i;
			TMPC_EXPECT_EQ(Y[i], Y_ref[i]) << "at i=" << i;
		}
	}


	TEST_F(GeneratedFunctionTest, testMatrixArgumentCall)
	{
		StaticMatrix<Real, 3, 2, columnMajor> A {
			{1., 2.},
			{3., 4.},
			{5., 6.}
		};
		StaticMatrix<Real, 2, 2, columnMajor> B {
			{7., 8.},
			{9., 10.}
		};
		Real x = 0.1;

		//randomize(A);
		//randomize(B);
		//blaze::randomize(x);

		StaticMatrix<Real, 3, 2> X;
		StaticVector<Real, 2, rowVector> Y;

		// fun_(std::tie(A, B, x));
		fun_(std::tie(A, B, x), std::tie(X, Y));

		TMPC_EXPECT_EQ(X, (A * x) * B);
		TMPC_EXPECT_EQ(Y, (StaticVector<Real, 3, rowVector> {1., 1., 1.} * (A * B)));
	}
}