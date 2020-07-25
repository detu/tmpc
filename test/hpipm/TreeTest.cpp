#include <tmpc/hpipm/Tree.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(TreeTest, testCtor)
	{
		OcpTree g {
			2,
			1, 1,
			0, 0
		};

		hpipm::Tree t {g};

		for (auto v : vertices(g))
		{
			ASSERT_EQ(t.root[v].nkids, out_degree(v, g));
			EXPECT_EQ(t.root[v].idx, v);
			EXPECT_EQ(t.root[v].stage, g.depth(v));
			EXPECT_EQ(t.root[v].dad, parent(v, g) ? *parent(v, g) : -1);

			for (int i = 0; i < t.root[v].nkids; ++i)
				EXPECT_EQ(t.root[v].kids[i], children(v, g)[i]);
		}
	}
}