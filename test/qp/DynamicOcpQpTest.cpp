/*
* Unit tests for common functionality of StaticOcpQp class.
*/

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    class DynamicOcpQpTest 
    :   public Test
    {
    public:
        using Real = double;

    protected:
        DynamicOcpQpTest()
        :   graph_(3)
        {
        }


        OcpTree const graph_;
        DynamicOcpSize const size_ {
            {3, 2, 1},
            {4, 3, 2},
            {3, 0, 1}
        };

        DynamicOcpQp<Real> qp_ {size_};
    };


    TEST_F(DynamicOcpQpTest, testQpInterface)
    {
        auto& g = this->qp_.graph();
        auto const N = num_vertices(g);

        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>>> Q;
        std::map<OcpVertex, blaze::DynamicVector<Real>> q;
        std::map<OcpVertex, blaze::DynamicVector<Real>> r;
        std::map<OcpVertex, blaze::DynamicMatrix<Real>> S;
        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>>> R;

        std::map<OcpEdge, blaze::DynamicMatrix<Real>> A;
        std::map<OcpEdge, blaze::DynamicMatrix<Real>> B;
        std::map<OcpEdge, blaze::DynamicVector<Real>> b;

        std::map<OcpVertex, blaze::DynamicVector<Real>> x_min, x_max;
        std::map<OcpVertex, blaze::DynamicVector<Real>> u_min, u_max;
        std::map<OcpVertex, blaze::DynamicVector<Real>> d_min, d_max;

        // Writing random data
        for (auto v : vertices(g))
        {
            Q[v].resize(size_.nx(v), size_.nx(v));
            R[v].resize(size_.nu(v), size_.nu(v));
            S[v].resize(size_.nu(v), size_.nx(v));
            q[v].resize(size_.nx(v));
            r[v].resize(size_.nu(v));
            x_min[v].resize(size_.nx(v));
            x_max[v].resize(size_.nx(v));
            u_min[v].resize(size_.nu(v));
            u_max[v].resize(size_.nu(v));
            d_min[v].resize(size_.nc(v));
            d_max[v].resize(size_.nc(v));

            makePositiveDefinite(Q[v]);
            randomize(q[v]);
            randomize(x_min[v]);
            randomize(x_max[v]);
            randomize(d_min[v]);
            randomize(d_max[v]);

            this->qp_.Q(v, Q[v]);
            this->qp_.q(v, q[v]);
            this->qp_.lx(v, x_min[v]);
            this->qp_.ux(v, x_max[v]);
            this->qp_.ld(v, d_min[v]);
            this->qp_.ud(v, d_max[v]);
        }

        for (auto v : g.branchVertices())
        {
            makePositiveDefinite(R[v]);
            randomize(S[v]);
            randomize(r[v]);
            randomize(u_min[v]);
            randomize(u_max[v]);            
        
            this->qp_.R(v, R[v]);
            this->qp_.S(v, S[v]);
            this->qp_.r(v, r[v]);
            this->qp_.lu(v, u_min[v]);
            this->qp_.uu(v, u_max[v]);
        }

        for (auto e : edges(g))
        {
            auto const u = source(e, graph_);
            auto const v = target(e, graph_);

            A[e].resize(size_.nx(v), size_.nx(u));
            B[e].resize(size_.nx(v), size_.nu(u));
            b[e].resize(size_.nx(v));

            randomize(A[e]);
            randomize(B[e]);
            randomize(b[e]);

            this->qp_.A(e, A[e]);
            this->qp_.B(e, B[e]);
            this->qp_.b(e, b[e]);
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : vertices(g))
        {
            EXPECT_EQ(forcePrint(this->qp_.Q(v)), forcePrint(Q[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.q(v)), forcePrint(q[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.lx(v)), forcePrint(x_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.ux(v)), forcePrint(x_max[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.ld(v)), forcePrint(d_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.ud(v)), forcePrint(d_max[v])) << "at v=" << v;
        }

        for (auto v : g.branchVertices())
        {
            EXPECT_EQ(forcePrint(this->qp_.R(v)), forcePrint(R[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.S(v)), forcePrint(S[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.r(v)), forcePrint(r[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.lu(v)), forcePrint(u_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(this->qp_.uu(v)), forcePrint(u_max[v])) << "at v=" << v;
        }

        for (auto e : edges(g))
        {
            EXPECT_EQ(forcePrint(this->qp_.A(e)), forcePrint(A[e])) << "at e=" << e;
            EXPECT_EQ(forcePrint(this->qp_.B(e)), forcePrint(B[e])) << "at e=" << e;
            EXPECT_EQ(forcePrint(this->qp_.b(e)), forcePrint(b[e])) << "at e=" << e;
        }
    }


    TEST_F(DynamicOcpQpTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->qp_;

        for (auto v : vertices(ws.graph()))
        {
            EXPECT_EQ(rows(this->qp_.Q(v)), size_.nx(v));
            EXPECT_EQ(columns(this->qp_.Q(v)), size_.nx(v));
            EXPECT_EQ(size(this->qp_.q(v)), size_.nx(v));
            EXPECT_EQ(size(this->qp_.lx(v)), size_.nx(v));
            EXPECT_EQ(size(this->qp_.ux(v)), size_.nx(v));
            EXPECT_EQ(size(this->qp_.ld(v)), size_.nc(v));
            EXPECT_EQ(size(this->qp_.ud(v)), size_.nc(v));
        }

        for (auto v : ws.graph().branchVertices())
        {
            EXPECT_EQ(rows(this->qp_.R(v)), size_.nu(v));
            EXPECT_EQ(columns(this->qp_.R(v)), size_.nu(v));
            EXPECT_EQ(rows(this->qp_.S(v)), size_.nu(v));
            EXPECT_EQ(columns(this->qp_.S(v)), size_.nx(v));
            EXPECT_EQ(size(this->qp_.r(v)), size_.nu(v));
            EXPECT_EQ(size(this->qp_.lu(v)), size_.nu(v));
            EXPECT_EQ(size(this->qp_.uu(v)), size_.nu(v));
        }

        for (auto e : edges(ws.graph()))
        {
            auto const u = source(e, graph_);
            auto const v = target(e, graph_);

            EXPECT_EQ(rows(this->qp_.A(e)), size_.nx(v));
            EXPECT_EQ(columns(this->qp_.A(e)), size_.nx(u));
            EXPECT_EQ(rows(this->qp_.B(e)), size_.nx(v));
            EXPECT_EQ(columns(this->qp_.B(e)), size_.nu(u));
            EXPECT_EQ(size(this->qp_.b(e)), size_.nx(v));
        }
    }
}
