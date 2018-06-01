/*
* Unit tests for common functionality of QP workspaces.
*/

#pragma once

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/GraphTools.hpp>
#include <tmpc/Matrix.hpp>

#include <tmpc/test_tools.hpp>
#include <gtest/gtest.h>

#include <iostream>
#include <array>

namespace tmpc :: testing
{
    template <typename WS>
    class TreeQpWorkspaceTest 
    :   public ::testing::Test
    {
    public:
        using Workspace = WS;

    protected:
        using Kernel = typename Workspace::Kernel;
        using Real = typename Workspace::Real;
        using Vector = DynamicVector<Kernel>;
        using Matrix = DynamicMatrix<Kernel>;


        TreeQpWorkspaceTest()
        :	size_ {initSizeGraph()}
        ,	ws_ {size_}
        {
        }


        OcpSizeGraph const size_;
        Workspace ws_;


        static OcpSizeGraph initSizeGraph()
        {
            OcpSizeGraph g(3);
            auto size = get(OcpSizeProperty_t(), g);

            put(size, 0, OcpSize(2, 3, 4));
            put(size, 1, OcpSize(5, 6, 7));
            put(size, 2, OcpSize(8, 9, 10));
            
            add_edge(0, 1, 0, g);
            add_edge(1, 2, 1, g);

            auto edge_index = get(boost::edge_index, g);
            for (auto e : edgesR(g))
            {
                auto const from = source(e, g);
                auto const to = target(e, g);
                std::cout << "from=" << from << ", to=" << to << ", edge_index=" << edge_index(e) << std::endl;
            }

            return g;
        }
    };


    TYPED_TEST_CASE_P(TreeQpWorkspaceTest);

    TYPED_TEST_P(TreeQpWorkspaceTest, testQpInterface)
    {
        /*
        auto const N = this->size_.size();

        std::vector<typename TestFixture::Matrix> Q(N);
        std::vector<typename TestFixture::Vector> q(N);
        std::vector<typename TestFixture::Vector> r(N);

        std::vector<typename TestFixture::Matrix> S(N);
        std::vector<typename TestFixture::Matrix> R(N);
        std::vector<typename TestFixture::Matrix> A(N);
        std::vector<typename TestFixture::Matrix> B(N);
        std::vector<typename TestFixture::Vector> b(N);

        std::vector<typename TestFixture::Vector> x_min(N), x_max(N);
        std::vector<typename TestFixture::Vector> u_min(N), u_max(N);
        std::vector<typename TestFixture::Vector> d_min(N), d_max(N);

        // Writing random data
        Rand<typename TestFixture::Kernel, typename TestFixture::Matrix> rand_matrix;
        Rand<typename TestFixture::Kernel, typename TestFixture::Vector> rand_vector;

        for (std::size_t i = 0; i < N; ++i)
        {
            auto const& sz = this->size_[i];
            auto const nx1 = i + 1 < N ? this->size_[i + 1].nx() : 0;
            auto& stage = this->ws_.problem()[i];

            stage.Q(Q[i] = rand_matrix.generate(sz.nx(), sz.nx()));
            stage.R(R[i] = rand_matrix.generate(sz.nu(), sz.nu()));
            stage.S(S[i] = rand_matrix.generate(sz.nx(), sz.nu()));
            stage.q(q[i] = rand_vector.generate(sz.nx()));
            stage.r(r[i] = rand_vector.generate(sz.nu()));

            stage.A(A[i] = rand_matrix.generate(nx1, sz.nx()));
            stage.B(B[i] = rand_matrix.generate(nx1, sz.nu()));
            stage.b(b[i] = rand_vector.generate(nx1));

            stage.lbx(x_min[i] = rand_vector.generate(sz.nx()));
            stage.ubx(x_max[i] = rand_vector.generate(sz.nx()));
            stage.lbu(u_min[i] = rand_vector.generate(sz.nu()));
            stage.ubu(u_max[i] = rand_vector.generate(sz.nu()));
            stage.lbd(d_min[i] = rand_vector.generate(sz.nc()));
            stage.ubd(d_max[i] = rand_vector.generate(sz.nc()));
        }

        // Reading the data and checking that they are the same that we wrote
        for (std::size_t i = 0; i < N; ++i)
        {
            auto const& stage = this->ws_.problem()[i];

            EXPECT_EQ(print_wrap(stage.Q()), print_wrap(Q[i])) << "at i=" << i;
            EXPECT_EQ(print_wrap(stage.R()), print_wrap(R[i])) << "at i=" << i;
            EXPECT_EQ(print_wrap(stage.S()), print_wrap(S[i]));
            EXPECT_EQ(print_wrap(stage.q()), print_wrap(q[i]));
            EXPECT_EQ(print_wrap(stage.r()), print_wrap(r[i]));

            EXPECT_EQ(print_wrap(stage.A()), print_wrap(A[i]));
            EXPECT_EQ(print_wrap(stage.B()), print_wrap(B[i]));
            EXPECT_EQ(print_wrap(stage.b()), print_wrap(b[i]));

            EXPECT_EQ(print_wrap(stage.lbx()), print_wrap(x_min[i]));
            EXPECT_EQ(print_wrap(stage.ubx()), print_wrap(x_max[i]));
            EXPECT_EQ(print_wrap(stage.lbu()), print_wrap(u_min[i]));
            EXPECT_EQ(print_wrap(stage.ubu()), print_wrap(u_max[i]));
            EXPECT_EQ(print_wrap(stage.lbd()), print_wrap(d_min[i]));
            EXPECT_EQ(print_wrap(stage.ubd()), print_wrap(d_max[i]));
        }
        */
    }

    TYPED_TEST_P(TreeQpWorkspaceTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->ws_;
        auto const original_size = get(OcpSizeProperty_t(), this->size_);

        for (auto v : verticesR(ws.graph()))
        {
            auto const& s = get(ws.size(), v);
            EXPECT_EQ(s, original_size(v));

            auto const& stage = get(ws.problemVertex(), v);

            EXPECT_EQ(rows   (stage.Q()), s.nx());
            EXPECT_EQ(columns(stage.Q()), s.nx());
            EXPECT_EQ(rows   (stage.R()), s.nu());
            EXPECT_EQ(columns(stage.R()), s.nu());
            EXPECT_EQ(rows   (stage.S()), s.nx());
            EXPECT_EQ(columns(stage.S()), s.nu());
            EXPECT_EQ(size   (stage.q()), s.nx());
            EXPECT_EQ(size   (stage.r()), s.nu());

            EXPECT_EQ(size(stage.lbx()), s.nx());
            EXPECT_EQ(size(stage.ubx()), s.nx());
            EXPECT_EQ(size(stage.lbu()), s.nu());
            EXPECT_EQ(size(stage.ubu()), s.nu());
            EXPECT_EQ(size(stage.lbd()), s.nc());
            EXPECT_EQ(size(stage.ubd()), s.nc());
        }

        auto edge_index = get(boost::edge_index, ws.graph());

        for (auto e : edgesR(ws.graph()))
        {
            auto const from = source(e, ws.graph());
            auto const to = target(e, ws.graph());

            auto const& sz_from = get(ws.size(), from);
            auto const& sz_to = get(ws.size(), to);

            std::cout << "from=" << from << ", to=" << to << ", edge_index=" << edge_index(e) << std::endl;

            auto const& problem_edge = get(ws.problemEdge(), e);

            EXPECT_EQ(rows(problem_edge.A()), sz_to.nx());
            EXPECT_EQ(columns(problem_edge.A()), sz_from.nx());
            EXPECT_EQ(rows(problem_edge.B()), sz_to.nx());
            EXPECT_EQ(columns(problem_edge.B()), sz_from.nu());
            EXPECT_EQ(size(problem_edge.b()), sz_to.nx());
        }
    }

    REGISTER_TYPED_TEST_CASE_P(TreeQpWorkspaceTest,
        testQpInterface, testMatrixSizesCorrect);
}