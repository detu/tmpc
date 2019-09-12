/*
* Unit tests for common functionality of QP workspaces.
*/

#pragma once

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/Traits.hpp>
#include <tmpc/Testing.hpp>

#include <iostream>
#include <array>


namespace tmpc :: testing
{
    template <typename WS>
    class TreeQpWorkspaceTest 
    :   public Test
    {
    public:
        using Workspace = WS;

    protected:
        using Kernel = typename KernelOf<Workspace>::type;
        using Real = typename RealOf<Workspace>::type;
        using Vector = DynamicVector<Kernel>;
        using Matrix = DynamicMatrix<Kernel>;


        TreeQpWorkspaceTest()
        {
        }


        OcpGraph const graph_ = ocpGraphLinear(3);

        std::vector<OcpSize> const size_ = {
            OcpSize {2, 1, 1},
            OcpSize {5, 4, 2},
            OcpSize {8, 7, 3}
        };

        Workspace ws_ {graph_, iterator_property_map(size_.begin(), vertexIndex(graph_))};
    };


    TYPED_TEST_SUITE_P(TreeQpWorkspaceTest);


    TYPED_TEST_P(TreeQpWorkspaceTest, testQpInterface)
    {
        auto& g = this->ws_.graph();
        auto const N = num_vertices(g);

        std::map<OcpVertexDescriptor, typename TestFixture::Matrix> Q;
        std::map<OcpVertexDescriptor, typename TestFixture::Vector> q;
        std::map<OcpVertexDescriptor, typename TestFixture::Vector> r;
        std::map<OcpVertexDescriptor, typename TestFixture::Matrix> S;
        std::map<OcpVertexDescriptor, typename TestFixture::Matrix> R;

        std::map<OcpEdgeDescriptor, typename TestFixture::Matrix> A;
        std::map<OcpEdgeDescriptor, typename TestFixture::Matrix> B;
        std::map<OcpEdgeDescriptor, typename TestFixture::Vector> b;

        std::map<OcpVertexDescriptor, typename TestFixture::Vector> x_min, x_max;
        std::map<OcpVertexDescriptor, typename TestFixture::Vector> u_min, u_max;
        std::map<OcpVertexDescriptor, typename TestFixture::Vector> d_min, d_max;

        // Writing random data
        Rand<typename TestFixture::Kernel, typename TestFixture::Matrix> rand_matrix;
        Rand<typename TestFixture::Kernel, typename TestFixture::Vector> rand_vector;

        for (auto v : graph::vertices(g))
        {
            auto const& sz = get(this->ws_.size(), v);

            Q[v].resize(sz.nx(), sz.nx());
            makePositiveDefinite(Q[v]);

            R[v].resize(sz.nu(), sz.nu());
            makePositiveDefinite(R[v]);

            put(this->ws_.Q(), v, Q[v]);
            put(this->ws_.R(), v, R[v]);
            put(this->ws_.S(), v, S[v] = rand_matrix.generate(sz.nu(), sz.nx()));
            put(this->ws_.q(), v, q[v] = rand_vector.generate(sz.nx()));
            put(this->ws_.r(), v, r[v] = rand_vector.generate(sz.nu()));

            put(this->ws_.lx(), v, x_min[v] = rand_vector.generate(sz.nx()));
            put(this->ws_.ux(), v, x_max[v] = rand_vector.generate(sz.nx()));
            put(this->ws_.lu(), v, u_min[v] = rand_vector.generate(sz.nu()));
            put(this->ws_.uu(), v, u_max[v] = rand_vector.generate(sz.nu()));
            put(this->ws_.ld(), v, d_min[v] = rand_vector.generate(sz.nc()));
            put(this->ws_.ud(), v, d_max[v] = rand_vector.generate(sz.nc()));
        }

        for (auto e : graph::edges(g))
        {
            auto const from = source(e, g);
            auto const to = target(e, g);
            auto const& sz_from = get(this->ws_.size(), from);
            auto const& sz_to = get(this->ws_.size(), to);

            put(this->ws_.A(), e, A[e] = rand_matrix.generate(sz_to.nx(), sz_from.nx()));
            put(this->ws_.B(), e, B[e] = rand_matrix.generate(sz_to.nx(), sz_from.nu()));
            put(this->ws_.b(), e, b[e] = rand_vector.generate(sz_to.nx()));
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : graph::vertices(g))
        {
            EXPECT_EQ(forcePrint(get(this->ws_.Q(), v)), forcePrint(Q[v])) << "at v=" << v;
            
            EXPECT_EQ(forcePrint(get(this->ws_.R(), v)), forcePrint(R[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.S(), v)), forcePrint(S[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.q(), v)), forcePrint(q[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.r(), v)), forcePrint(r[v])) << "at v=" << v;

            EXPECT_EQ(forcePrint(get(this->ws_.lx(), v)), forcePrint(x_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.ux(), v)), forcePrint(x_max[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.lu(), v)), forcePrint(u_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.uu(), v)), forcePrint(u_max[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.ld(), v)), forcePrint(d_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->ws_.ud(), v)), forcePrint(d_max[v])) << "at v=" << v;
        }

        for (auto e : graph::edges(g))
        {
            EXPECT_EQ(forcePrint(get(this->ws_.A(), e)), forcePrint(A[e])) << "at e=" << get(graph::edge_index, g, e);
            EXPECT_EQ(forcePrint(get(this->ws_.B(), e)), forcePrint(B[e])) << "at e=" << get(graph::edge_index, g, e);
            EXPECT_EQ(forcePrint(get(this->ws_.b(), e)), forcePrint(b[e])) << "at e=" << get(graph::edge_index, g, e);
        }
    }


    TYPED_TEST_P(TreeQpWorkspaceTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->ws_;

        for (auto v : graph::vertices(ws.graph()))
        {
            auto const& s = get(ws.size(), v);
            EXPECT_EQ(s, this->size_[v]);

            EXPECT_EQ(rows(get(this->ws_.Q(), v)), s.nx());
            EXPECT_EQ(columns(get(this->ws_.Q(), v)), s.nx());

            EXPECT_EQ(rows(get(this->ws_.R(), v)), s.nu());
            EXPECT_EQ(columns(get(this->ws_.R(), v)), s.nu());
            EXPECT_EQ(rows(get(this->ws_.S(), v)), s.nu());
            EXPECT_EQ(columns(get(this->ws_.S(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.q(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.r(), v)), s.nu());

            EXPECT_EQ(size(get(this->ws_.lx(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.ux(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.lu(), v)), s.nu());
            EXPECT_EQ(size(get(this->ws_.uu(), v)), s.nu());
            EXPECT_EQ(size(get(this->ws_.ld(), v)), s.nc());
            EXPECT_EQ(size(get(this->ws_.ud(), v)), s.nc());
        }

        for (auto e : graph::edges(ws.graph()))
        {
            auto const from = source(e, ws.graph());
            auto const to = target(e, ws.graph());

            auto const& sz_from = get(ws.size(), from);
            auto const& sz_to = get(ws.size(), to);

            //std::cout << "from=" << from << ", to=" << to << ", edge_index=" << edge_index(e) << std::endl;

            //auto const& problem_edge = get(ws.problemEdge(), e);

            EXPECT_EQ(rows(get(this->ws_.A(), e)), sz_to.nx());
            EXPECT_EQ(columns(get(this->ws_.A(), e)), sz_from.nx());
            EXPECT_EQ(rows(get(this->ws_.B(), e)), sz_to.nx());
            EXPECT_EQ(columns(get(this->ws_.B(), e)), sz_from.nu());
            EXPECT_EQ(size(get(this->ws_.b(), e)), sz_to.nx());
        }
    }

    REGISTER_TYPED_TEST_SUITE_P(TreeQpWorkspaceTest,
        testQpInterface, testMatrixSizesCorrect);
}
