#include <tmpc/casadi_interface/GeneratedFunction.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>
#include <tmpc/Testing.hpp>

#include <test_functions.h>


#include <tuple>
#include <stdexcept>
#include <random>


namespace tmpc :: testing
{
	using Kernel = EigenKernel<double>;
	using namespace casadi_interface;


	TEST(CompressedColumnStorageToMatrixTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		compressedColumnStorageToMatrix(data, sparsity, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		}));
	}


	TEST(CompressedColumnStorageToMatrixTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real const data[3] = {1.1, 2.2, 3.3};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		compressedColumnStorageToMatrix(data, sparsity, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		}));
	}


	TEST(CompressedColumnStorageToMatrixTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> m;
		EXPECT_THROW(compressedColumnStorageToMatrix(data, sparsity, m), std::invalid_argument);
	}


	TEST(CompressedColumnStorageToMatrixTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> m;
		EXPECT_THROW(compressedColumnStorageToMatrix(data, sparsity, m), std::invalid_argument);
	}


	TEST(CompressedColumnStorageToMatrixTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		EXPECT_THROW(compressedColumnStorageToMatrix(data, sparsity, m), std::invalid_argument);
	}


	TEST(MatrixToCompressedColumnStorageTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {0., 0., 0., 0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		};

		matrixToCompressedColumnStorage(m, data, sparsity);

		EXPECT_THAT(data, ElementsAre(1.1, 2.1, 3.1, 1.2, 2.2, 3.2));
	}


	TEST(MatrixToCompressedColumnStorageTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real data[3] = {0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		};

		matrixToCompressedColumnStorage(m, data, sparsity);

		EXPECT_THAT(data, ElementsAre(1.1, 2.2, 3.3));
	}


	TEST(MatrixToCompressedColumnStorageTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> const m(0.);
		EXPECT_THROW(matrixToCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}


	TEST(MatrixToCompressedColumnStorageTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> const m(0.);
		EXPECT_THROW(matrixToCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}


	TEST(MatrixToCompressedColumnStorageTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> const m(0.);
		EXPECT_THROW(matrixToCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}


	class GeneratedFunctionTest 
	: 	public Test
	{
	protected:
		casadi_interface::GeneratedFunction fun_ {f_functions()};
		using MatrixA = StaticMatrix<Kernel, 3, 2, columnMajor>;
		using MatrixB = StaticMatrix<Kernel, 2, 2, columnMajor>;
		using MatrixX = StaticMatrix<Kernel, 3, 2>;
		using VectorY = StaticVector<Kernel, 2, rowVector>;
		using Real = Kernel::Real;
	};


	TEST_F(GeneratedFunctionTest, incorrect_n_inputs_throws)
	{
		Kernel::Real x = 0.;
		ASSERT_THROW(fun_({&x, &x, &x, &x}, {&x, &x}), std::invalid_argument);
	}


	TEST_F(GeneratedFunctionTest, incorrect_n_outputs_throws)
	{
		Kernel::Real x = 0.;
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
		StaticMatrix<Kernel, 3, 2, columnMajor> A {
			{1., 2.},
			{3., 4.},
			{5., 6.}
		};
		StaticMatrix<Kernel, 2, 2, columnMajor> B {
			{7., 8.},
			{9., 10.}
		};
		Kernel::Real x = 0.1;

		//randomize(A);
		//randomize(B);
		//blaze::randomize(x);

		StaticMatrix<Kernel, 3, 2> X;
		StaticVector<Kernel, 2, rowVector> Y;

		fun_({A.data(), B.data(), &x}, {X.data(), Y.data()});

		EXPECT_PRED2(ApproxEqual(1e-12), X, (A * x) * B);
		EXPECT_EQ(Y, (StaticVector<Kernel, 3, rowVector> {1., 1., 1.} * (A * B)));
	}


	TEST_F(GeneratedFunctionTest, testCopyCtor)
	{
		StaticMatrix<Kernel, 3, 2, columnMajor> A {
			{1., 2.},
			{3., 4.},
			{5., 6.}
		};
		StaticMatrix<Kernel, 2, 2, columnMajor> B {
			{7., 8.},
			{9., 10.}
		};
		Kernel::Real x = 0.1;

		//randomize(A);
		//randomize(B);
		//blaze::randomize(x);

		StaticMatrix<Kernel, 3, 2> X, X1;
		StaticVector<Kernel, 2, rowVector> Y, Y1;

		casadi_interface::GeneratedFunction const fun_copy = fun_;
		fun_({A.data(), B.data(), &x}, {X.data(), Y.data()});
		fun_copy({A.data(), B.data(), &x}, {X1.data(), Y1.data()});

		EXPECT_PRED2(ApproxEqual(1e-12), X, (A * x) * B);
		EXPECT_EQ(Y, (StaticVector<Kernel, 3, rowVector> {1., 1., 1.} * (A * B)));
		EXPECT_PRED2(ApproxEqual(1e-12), X1, (A * x) * B);
		EXPECT_EQ(Y1, (StaticVector<Kernel, 3, rowVector> {1., 1., 1.} * (A * B)));
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
			A[i] = MatrixA::Random();
			B[i] = MatrixB::Random();
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
			EXPECT_EQ(forcePrint(X[i]), forcePrint(X_ref[i])) << "at i=" << i;
			EXPECT_EQ(forcePrint(Y[i]), forcePrint(Y_ref[i])) << "at i=" << i;
		}
	}

	/*
	TEST_F(GeneratedFunctionTest, testMatrixArgumentCall)
	{
		StaticMatrix<Kernel, 3, 2, columnMajor> A {
			{1., 2.},
			{3., 4.},
			{5., 6.}
		};
		StaticMatrix<Kernel, 2, 2, columnMajor> B {
			{7., 8.},
			{9., 10.}
		};
		Kernel::Real x = 0.1;

		//randomize(A);
		//randomize(B);
		//blaze::randomize(x);

		StaticMatrix<Kernel, 3, 2> X;
		StaticVector<Kernel, 2, rowVector> Y;

		fun_({A, B, x}, {X, Y});

		EXPECT_EQ(forcePrint(X), forcePrint((A * x) * B));
		EXPECT_EQ(Y, (StaticVector<Kernel, 3, rowVector> {1., 1., 1.} * (A * B)));
	}
	*/

	/*
	// Some interesting ideas here: https://habrahabr.ru/post/228031/
	template <typename Function, typename TupleIn, typename TupleOut>
	void call(Function const& f, TupleIn&& in, TupleOut&& out)
	{
		throw std::logic_error("Not implemented");
	}

	TEST_F(GeneratedFunctionTest, DISABLED_matrix_argument_call_correct)
	{
		StaticMatrix<Kernel, 3, 2> const A {0.};
		StaticMatrix<Kernel, 2, 2> const B {0.};
		StaticMatrix<Kernel, 1, 1> const x {0.};

		StaticMatrix<Kernel, 3, 2> X;
		StaticMatrix<Kernel, 1, 2> Y;

		call(fun_, std::forward_as_tuple(A, B, x), std::forward_as_tuple(X, Y));
	}
	*/
}