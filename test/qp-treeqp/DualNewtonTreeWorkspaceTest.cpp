#include "TreeQpWorkspaceSolveTest.hpp"

//#include "TreeQpWorkspaceTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>

/*
* Unit tests for common functionality of QP workspaces.
*/

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

            g[0].size = OcpSize(2, 1, 1);
            g[1].size = OcpSize(5, 4, 2);
            g[2].size = OcpSize(8, 7, 3);
            
            add_edge(0, 1, 0, g);
            add_edge(1, 2, 1, g);

            return g;
        }
    };


    TYPED_TEST_CASE_P(TreeQpWorkspaceTest);


    TYPED_TEST_P(TreeQpWorkspaceTest, testQpInterface)
    {
        using graph_traits = boost::graph_traits<OcpSizeGraph>;

        auto const N = num_vertices(this->ws_.graph());

        std::map<graph_traits::vertex_descriptor, typename TestFixture::Matrix> Q;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Vector> q;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Vector> r;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Matrix> S;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Matrix> R;

        std::map<graph_traits::edge_descriptor, typename TestFixture::Matrix> A;
        std::map<graph_traits::edge_descriptor, typename TestFixture::Matrix> B;
        std::map<graph_traits::edge_descriptor, typename TestFixture::Vector> b;

        std::map<graph_traits::vertex_descriptor, typename TestFixture::Vector> x_min, x_max;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Vector> u_min, u_max;
        std::map<graph_traits::vertex_descriptor, typename TestFixture::Vector> d_min, d_max;

        // Writing random data
        Rand<typename TestFixture::Kernel, typename TestFixture::Matrix> rand_matrix;
        Rand<typename TestFixture::Kernel, typename TestFixture::Vector> rand_vector;

        for (auto v : verticesR(this->ws_.graph()))
        {
            auto const& sz = get(this->ws_.size(), v);
            //auto& stage = get(this->ws_.problemVertex(), v);

            put(this->ws_.Q(), v, Q[v] = rand_matrix.generate(sz.nx(), sz.nx()));
            /*
            stage.R(R[v] = rand_matrix.generate(sz.nu(), sz.nu()));
            stage.S(S[v] = rand_matrix.generate(sz.nx(), sz.nu()));
            stage.q(q[v] = rand_vector.generate(sz.nx()));
            stage.r(r[v] = rand_vector.generate(sz.nu()));

            stage.lbx(x_min[v] = rand_vector.generate(sz.nx()));
            stage.ubx(x_max[v] = rand_vector.generate(sz.nx()));
            stage.lbu(u_min[v] = rand_vector.generate(sz.nu()));
            stage.ubu(u_max[v] = rand_vector.generate(sz.nu()));
            stage.lbd(d_min[v] = rand_vector.generate(sz.nc()));
            stage.ubd(d_max[v] = rand_vector.generate(sz.nc()));
            */
        }

        for (auto e : edgesR(this->ws_.graph()))
        {
            auto const from = source(e, this->ws_.graph());
            auto const to = target(e, this->ws_.graph());
            auto const& sz_from = get(this->ws_.size(), from);
            auto const& sz_to = get(this->ws_.size(), to);

            /*
            auto& stage = get(this->ws_.problemEdge(), e);

            stage.A(A[e] = rand_matrix.generate(sz_to.nx(), sz_from.nx()));
            stage.B(B[e] = rand_matrix.generate(sz_to.nx(), sz_from.nu()));
            stage.b(b[e] = rand_vector.generate(sz_to.nx()));
            */
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : verticesR(this->ws_.graph()))
        {
            //auto const& stage = get(this->ws_.problemVertex(), v);

            EXPECT_EQ(print_wrap(get(this->ws_.Q(), v)), print_wrap(Q[v])) << "at v=" << v;
            
            /*
            EXPECT_EQ(print_wrap(stage.R()), print_wrap(R[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.S()), print_wrap(S[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.q()), print_wrap(q[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.r()), print_wrap(r[v])) << "at v=" << v;

            EXPECT_EQ(print_wrap(stage.lbx()), print_wrap(x_min[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.ubx()), print_wrap(x_max[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.lbu()), print_wrap(u_min[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.ubu()), print_wrap(u_max[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.lbd()), print_wrap(d_min[v])) << "at v=" << v;
            EXPECT_EQ(print_wrap(stage.ubd()), print_wrap(d_max[v])) << "at v=" << v;
            */
        }

        for (auto e : edgesR(this->ws_.graph()))
        {
            /*
            auto const& stage = get(this->ws_.problemEdge(), e);
            EXPECT_EQ(print_wrap(stage.A()), print_wrap(A[e])) << "at e=" << e;
            EXPECT_EQ(print_wrap(stage.B()), print_wrap(B[e])) << "at e=" << e;
            EXPECT_EQ(print_wrap(stage.b()), print_wrap(b[e])) << "at e=" << e;
            */
        }
    }


    TYPED_TEST_P(TreeQpWorkspaceTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->ws_;

        for (auto v : verticesR(ws.graph()))
        {
            auto const& s = get(ws.size(), v);
            EXPECT_EQ(s, this->size_[v].size);

            EXPECT_EQ(rows(get(this->ws_.Q(), v)), s.nx());
            EXPECT_EQ(columns(get(this->ws_.Q(), v)), s.nx());

            EXPECT_EQ(rows(get(this->ws_.R(), v)), s.nu());
            EXPECT_EQ(columns(get(this->ws_.R(), v)), s.nu());
            EXPECT_EQ(rows(get(this->ws_.S(), v)), s.nx());
            EXPECT_EQ(columns(get(this->ws_.S(), v)), s.nu());
            EXPECT_EQ(size(get(this->ws_.q(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.r(), v)), s.nu());

            EXPECT_EQ(size(get(this->ws_.lx(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.ux(), v)), s.nx());
            EXPECT_EQ(size(get(this->ws_.lu(), v)), s.nu());
            EXPECT_EQ(size(get(this->ws_.uu(), v)), s.nu());
            EXPECT_EQ(size(get(this->ws_.ld(), v)), s.nc());
            EXPECT_EQ(size(get(this->ws_.ud(), v)), s.nc());
        }

        auto edge_index = get(boost::edge_index, ws.graph());

        for (auto e : edgesR(ws.graph()))
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

    REGISTER_TYPED_TEST_CASE_P(TreeQpWorkspaceTest,
        testQpInterface, testMatrixSizesCorrect);
}


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Blaze_double, TreeQpWorkspaceTest, DualNewtonTreeWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Eigen_double, TreeQpWorkspaceTest, DualNewtonTreeWorkspace<EigenKernel<double>>);
    INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Blaze_double, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace<BlazeKernel<double>>);
    INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Eigen_double, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace<EigenKernel<double>>);
}