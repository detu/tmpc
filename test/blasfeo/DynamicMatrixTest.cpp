#include <tmpc/blasfeo/DynamicMatrix.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(DynamicMatrixTest, testCtor)
	{
		blasfeo::DynamicMatrix<double> m(2, 3);
	}


	TEST(DynamicMatrixTest, testRows)
	{
		blasfeo::DynamicMatrix<double> m(2, 3);
		EXPECT_EQ(rows(m), 2);
	}


	TEST(DynamicMatrixTest, testColumns)
	{
		blasfeo::DynamicMatrix<double> m(2, 3);
		EXPECT_EQ(columns(m), 3);
	}


	TEST(DynamicMatrixTest, testResize)
	{
		blasfeo::DynamicMatrix<double> m(2, 3);
		ASSERT_EQ(rows(m), 2);
		ASSERT_EQ(columns(m), 3);

		m.resize(4, 5);
		ASSERT_EQ(rows(m), 4);
		ASSERT_EQ(columns(m), 5);

		// Assign to the last element to check that there is enough memory s.t. we don't SEGFAULT
		m(3, 4) = 42.;
	}


	TEST(DynamicMatrixTest, testElementAccess)
	{
		blasfeo::DynamicMatrix<double> m(2, 3);
		ASSERT_EQ(rows(m), 2);
		ASSERT_EQ(columns(m), 3);

		for (size_t i = 0; i < rows(m); ++i)
			for (size_t j = 0; j < columns(m); ++j)
				m(i, j) = 10 * i + j;

		for (size_t i = 0; i < rows(m); ++i)
			for (size_t j = 0; j < columns(m); ++j)
				EXPECT_EQ(m(i, j), 10 * i + j) << "element mismatch at index (" << i << ", " << j << ")";
	}
}