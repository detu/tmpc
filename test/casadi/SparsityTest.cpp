#include <tmpc/casadi/Sparsity.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>

#include <stdexcept>


namespace tmpc :: testing
{
	using namespace tmpc :: casadi;


	TEST(CcsDecompressTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		Sparsity(sparsity).decompress(data, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		}));
	}


	TEST(CcsDecompressTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real const data[3] = {1.1, 2.2, 3.3};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		Sparsity(sparsity).decompress(data, m);

		EXPECT_EQ(forcePrint(m), forcePrint(blaze::StaticMatrix<casadi_real, 3, 2> {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		}));
	}


	TEST(CcsDecompressTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> m;
		EXPECT_THROW(Sparsity(sparsity).decompress(data, m), std::invalid_argument);
	}


	TEST(CcsDecompressTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> m;
		EXPECT_THROW(Sparsity(sparsity).decompress(data, m), std::invalid_argument);
	}


	TEST(CcsDecompressTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real const data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> m;
		EXPECT_THROW(Sparsity(sparsity).decompress(data, m), std::invalid_argument);
	}


	TEST(CcsCompressTest, testDense)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {0., 0., 0., 0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{1.1, 1.2},
			{2.1, 2.2},
			{3.1, 3.2}
		};

		Sparsity(sparsity).compress(m, data);

		EXPECT_THAT(data, ElementsAre(1.1, 2.1, 3.1, 1.2, 2.2, 3.2));
	}


	TEST(CcsCompressTest, testSparse)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 2, 3, 1, 2, 0};
		casadi_real data[3] = {0., 0., 0.};

		blaze::StaticMatrix<casadi_real, 3, 2> const m {
			{0.0, 3.3},
			{1.1, 0.0},
			{2.2, 0.0}
		};

		Sparsity(sparsity).compress(m, data);

		EXPECT_THAT(data, ElementsAre(1.1, 2.2, 3.3));
	}


	TEST(CcsCompressTest, testInvalidNumRowsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 4, 2> const m(0.);
		EXPECT_THROW(Sparsity(sparsity).compress(m, data), std::invalid_argument);
	}


	TEST(CcsCompressTest, testInvalidNumColsThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 0, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 3> const m(0.);
		EXPECT_THROW(Sparsity(sparsity).compress(m, data), std::invalid_argument);
	}


	TEST(CcsCompressTest, testInvalidSparsityThrows)
	{
		casadi_int const sparsity[11] = {3, 2, 1, 3, 6, 0, 1, 2, 0, 1, 2};
		casadi_real data[6] = {1.1, 2.1, 3.1, 1.2, 2.2, 3.2};

		blaze::StaticMatrix<casadi_real, 3, 2> const m(0.);
		EXPECT_THROW(Sparsity(sparsity).compress(m, data), std::invalid_argument);
	}
}