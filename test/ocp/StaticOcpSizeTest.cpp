#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(StaticOcpSizeTest, testSize)
	{
		static size_t constexpr NX = 3;
		static size_t constexpr NU = 2;
		static size_t constexpr NC = 1;

		StaticOcpSize<NX, NU, NC> const s {OcpTree {2, 1, 1, 1, 1, 0, 0}};

		EXPECT_EQ(s.nx(0), 3);
		EXPECT_EQ(s.nu(0), 2);
		EXPECT_EQ(s.nc(0), 1);
		
		EXPECT_EQ(s.nx(1), 3);
		EXPECT_EQ(s.nu(1), 2);
		EXPECT_EQ(s.nc(1), 1);

		EXPECT_EQ(s.nx(2), 3);
		EXPECT_EQ(s.nu(2), 2);
		EXPECT_EQ(s.nc(2), 1);

		EXPECT_EQ(s.nx(3), 3);
		EXPECT_EQ(s.nu(3), 2);
		EXPECT_EQ(s.nc(3), 1);

		EXPECT_EQ(s.nx(4), 3);
		EXPECT_EQ(s.nu(4), 2);
		EXPECT_EQ(s.nc(4), 1);

		EXPECT_EQ(s.nx(5), 3);
		EXPECT_EQ(s.nu(5), 0);
		EXPECT_EQ(s.nc(5), 1);

		EXPECT_EQ(s.nx(6), 3);
		EXPECT_EQ(s.nu(6), 0);
		EXPECT_EQ(s.nc(6), 1);
	}


	TEST(StaticOcpSizeTest, testEqual)
	{
		static size_t constexpr NX = 3;
		static size_t constexpr NU = 2;
		static size_t constexpr NC = 1;

		StaticOcpSize<NX, NU, NC> const s1 {OcpTree {2, 1, 1, 1, 1, 0, 0}};
		StaticOcpSize<NX, NU, NC> const s2 {OcpTree {2, 1, 1, 1, 1, 0, 0}};
		StaticOcpSize<NX, NU, NC> const s3 {OcpTree {2, 1, 1, 0, 0}};
		StaticOcpSize<NX + 1, NU, NC> const s4 {OcpTree {2, 1, 1, 1, 1, 0, 0}};

		EXPECT_EQ(s1, s1);
		EXPECT_EQ(s1, s2);
		EXPECT_EQ(s2, s1);
		EXPECT_NE(s1, s3);
		EXPECT_NE(s3, s1);
		EXPECT_NE(s4, s1);
	}
}