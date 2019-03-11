#include <tmpc/casadi_interface/CompressedColumnStorage.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>

#include <stdexcept>


namespace tmpc :: testing
{
	using namespace casadi_interface;


	TEST(FromCompressedColumnStorageTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		fromCompressedColumnStorage(data, sparsity, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		}));
	}


	TEST(FromCompressedColumnStorageTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real const data[3] = {1.1, 2.2, 3.3};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		fromCompressedColumnStorage(data, sparsity, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		}));
	}


	TEST(FromCompressedColumnStorageTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> m;
		EXPECT_THROW(fromCompressedColumnStorage(data, sparsity, m), std::invalid_argument);
	}


	TEST(FromCompressedColumnStorageTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> m;
		EXPECT_THROW(fromCompressedColumnStorage(data, sparsity, m), std::invalid_argument);
	}


	TEST(FromCompressedColumnStorageTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		EXPECT_THROW(fromCompressedColumnStorage(data, sparsity, m), std::invalid_argument);
	}


	TEST(ToCompressedColumnStorageTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {0., 0., 0., 0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		};

		toCompressedColumnStorage(m, data, sparsity);

		EXPECT_THAT(data, ElementsAre(1.1, 2.1, 3.1, 1.2, 2.2, 3.2));
	}


	TEST(ToCompressedColumnStorageTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real data[3] = {0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		};

		toCompressedColumnStorage(m, data, sparsity);

		EXPECT_THAT(data, ElementsAre(1.1, 2.2, 3.3));
	}


	TEST(ToCompressedColumnStorageTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> const m(0.);
		EXPECT_THROW(toCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}


	TEST(ToCompressedColumnStorageTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> const m(0.);
		EXPECT_THROW(toCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}


	TEST(ToCompressedColumnStorageTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> const m(0.);
		EXPECT_THROW(toCompressedColumnStorage(m, data, sparsity), std::invalid_argument);
	}
}