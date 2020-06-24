#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(CondensedOcpSizeTest, testInitializerListArg)
	{
		EXPECT_EQ(condensedOcpSize({
			OcpSize(2, 1, 3),
			OcpSize(4, 5, 6),
			OcpSize(2, 1, 1)
			}),
			OcpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
	}

	TEST(CondensedOcpSizeTest, testIteratorRangeArg)
	{
		std::array<OcpSize, 3> sz = {
			OcpSize(2, 1, 3),
			OcpSize(4, 5, 6),
			OcpSize(2, 1, 1)
			};

		EXPECT_EQ(condensedOcpSize(sz.begin(), sz.end()),
			OcpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
	}
}