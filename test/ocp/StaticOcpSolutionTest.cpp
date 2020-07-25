/*
* Unit tests for common functionality of StaticOcpQp class.
*/

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/ocp/StaticOcpSolution.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    class StaticOcpSolutionTest 
    :   public Test
    {
    public:
        using Real = double;
        static size_t constexpr NX = 4;
        static size_t constexpr NU = 3;
        static size_t constexpr NC = 2;

        using Solution = StaticOcpSolution<Real, NX, NU, NC>;

    protected:
        StaticOcpSolutionTest()
        :   graph_(3)
        {
        }


        OcpTree const graph_;
        Solution sol_ {graph_};
    };


    TEST_F(StaticOcpSolutionTest, testQpSolutionInterface)
    {
        std::map<OcpVertex, blaze::StaticVector<Real, NX>> x, lam_lx, lam_ux;
        std::map<OcpVertex, blaze::StaticVector<Real, NU>> u, lam_lu, lam_uu;
        std::map<OcpVertex, blaze::StaticVector<Real, NC>> lam_ld, lam_ud;

        std::map<OcpEdge, blaze::StaticVector<Real, NX>> pi;

        // Writing random data
        for (auto v : vertices(graph_))
        {
            randomize(x[v]);
            randomize(lam_lx[v]);
            randomize(lam_ux[v]);
            randomize(lam_ld[v]);
            randomize(lam_ud[v]);

            sol_.x(v) = x[v];
            sol_.lam_lx(v) = lam_lx[v];
            sol_.lam_ux(v) = lam_ux[v];
            sol_.lam_ld(v) = lam_ld[v];
            sol_.lam_ud(v) = lam_ud[v];
        }

        for (auto v : graph_.branchVertices())
        {
            randomize(u[v]);
            randomize(lam_lu[v]);
            randomize(lam_uu[v]);
            
            sol_.u(v) = u[v];
            sol_.lam_lu(v) = lam_lu[v];
            sol_.lam_uu(v) = lam_uu[v];
        }

        for (auto e : edges(graph_))
        {
            randomize(pi[e]);            
            sol_.pi(e) = pi[e];
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : vertices(graph_))
        {
            EXPECT_EQ(forcePrint(sol_.x(v)), forcePrint(x[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_lx(v)), forcePrint(lam_lx[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_ux(v)), forcePrint(lam_ux[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_ld(v)), forcePrint(lam_ld[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_ud(v)), forcePrint(lam_ud[v])) << "at v=" << v;
        }

        for (auto v : graph_.branchVertices())
        {
            EXPECT_EQ(forcePrint(sol_.u(v)), forcePrint(u[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_lu(v)), forcePrint(lam_lu[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(sol_.lam_uu(v)), forcePrint(lam_uu[v])) << "at v=" << v;
        }

        for (auto e : edges(graph_))
        {
            EXPECT_EQ(forcePrint(sol_.pi(e)), forcePrint(pi[e])) << "at e=" << e;
        }
    }


    TEST_F(StaticOcpSolutionTest, testMatrixSizesCorrect)
    {
        for (auto v : vertices(graph_))
        {
            EXPECT_EQ(size(sol_.x(v)), NX);
            EXPECT_EQ(size(sol_.lam_lx(v)), NX);
            EXPECT_EQ(size(sol_.lam_ux(v)), NX);

            EXPECT_EQ(size(sol_.lam_ld(v)), NC);
            EXPECT_EQ(size(sol_.lam_ud(v)), NC);
        }

        for (auto v : graph_.branchVertices())
        {
            EXPECT_EQ(size(sol_.u(v)), NU);
            EXPECT_EQ(size(sol_.lam_lu(v)), NU);
            EXPECT_EQ(size(sol_.lam_uu(v)), NU);
        }

        for (auto e : edges(graph_))
        {
            EXPECT_EQ(size(sol_.pi(e)), NX);
        }
    }
}
