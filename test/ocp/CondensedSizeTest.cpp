#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(CondensedSizeTest, testCondensedSize)
	{
		EXPECT_EQ(condensedSize(DynamicOcpSize {
			{2, 1, 3},
			{4, 5, 6},
			{2, 1, 1}
			}),
			(DynamicOcpSize {{2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)}})
		);
	}
}