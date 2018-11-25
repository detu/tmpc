#include <tmpc/blasfeo/blasfeo.hpp>

#include <tmpc/test_tools.hpp>

#include <vector>


namespace tmpc :: testing
{
	TEST(BlasfeoCustomMatrixTest, testCtor)
	{
		std::vector<double> data(6);
		blasfeo::CustomMatrix<double> m(data.data(), 2, 3);
	}


	TEST(BlasfeoCustomMatrixTest, testRows)
	{
		std::vector<double> data(6);
		blasfeo::CustomMatrix<double> m(data.data(), 2, 3);

		EXPECT_EQ(rows(m), 2);
	}


	TEST(BlasfeoCustomMatrixTest, testColumns)
	{
		std::vector<double> data(6);
		blasfeo::CustomMatrix<double> m(data.data(), 2, 3);

		EXPECT_EQ(columns(m), 3);
	}
}