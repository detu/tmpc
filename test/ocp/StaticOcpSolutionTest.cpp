/*
* Unit tests for common functionality of StaticOcpQp class.
*/

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/PropertyMap.hpp>
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
        {
        }


        OcpGraph const graph_ = ocpGraphLinear(3);
        Solution sol_ {graph_};
    };


    TEST_F(StaticOcpSolutionTest, testQpSolutionInterface)
    {
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NX>> x, lam_lx, lam_ux;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NU>> u, lam_lu, lam_uu;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NC>> lam_ld, lam_ud;

        std::map<OcpEdgeDescriptor, blaze::StaticVector<Real, NX>> pi;

        // Writing random data
        for (auto v : graph::vertices(graph_))
        {
            randomize(x[v]);
            randomize(lam_lx[v]);
            randomize(lam_ux[v]);
            randomize(u[v]);
            randomize(lam_lu[v]);
            randomize(lam_uu[v]);
            randomize(lam_ld[v]);
            randomize(lam_ud[v]);

            put(sol_.x(), v, x[v]);
            put(sol_.lam_lx(), v, lam_lx[v]);
            put(sol_.lam_ux(), v, lam_ux[v]);
            put(sol_.u(), v, u[v]);
            put(sol_.lam_lu(), v, lam_lu[v]);
            put(sol_.lam_uu(), v, lam_uu[v]);
            put(sol_.lam_ld(), v, lam_ld[v]);
            put(sol_.lam_ud(), v, lam_ud[v]);
        }

        for (auto e : graph::edges(graph_))
        {
            randomize(pi[e]);
            
            put(sol_.pi(), e, pi[e]);
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : graph::vertices(graph_))
        {
            EXPECT_EQ(forcePrint(get(sol_.x(), v)), forcePrint(x[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(sol_.lam_lx(), v)), forcePrint(lam_lx[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(sol_.lam_ux(), v)), forcePrint(lam_ux[v])) << "at v=" << v;

            EXPECT_EQ(forcePrint(get(sol_.u(), v)), forcePrint(u[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(sol_.lam_lu(), v)), forcePrint(lam_lu[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(sol_.lam_uu(), v)), forcePrint(lam_uu[v])) << "at v=" << v;

            EXPECT_EQ(forcePrint(get(sol_.lam_ld(), v)), forcePrint(lam_ld[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(sol_.lam_ud(), v)), forcePrint(lam_ud[v])) << "at v=" << v;
        }

        for (auto e : graph::edges(graph_))
        {
            EXPECT_EQ(forcePrint(get(sol_.pi(), e)), forcePrint(pi[e])) << "at e=" << get(graph::edge_index, graph_, e);
        }
    }


    TEST_F(StaticOcpSolutionTest, testMatrixSizesCorrect)
    {
        for (auto v : graph::vertices(graph_))
        {
            OcpSize const s {NX, NU, NC};
            
            EXPECT_EQ(size(get(sol_.x(), v)), s.nx());
            EXPECT_EQ(size(get(sol_.lam_lx(), v)), s.nx());
            EXPECT_EQ(size(get(sol_.lam_ux(), v)), s.nx());

            EXPECT_EQ(size(get(sol_.u(), v)), s.nu());
            EXPECT_EQ(size(get(sol_.lam_lu(), v)), s.nu());
            EXPECT_EQ(size(get(sol_.lam_uu(), v)), s.nu());

            EXPECT_EQ(size(get(sol_.lam_ld(), v)), s.nc());
            EXPECT_EQ(size(get(sol_.lam_ud(), v)), s.nc());
        }

        for (auto e : graph::edges(graph_))
        {
            OcpSize const s {NX, NU, NC};

            EXPECT_EQ(size(get(sol_.pi(), e)), s.nx());
        }
    }
}
