/*
* Unit tests for common functionality of StaticOcpQp class.
*/

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/qp/StaticOcpQp.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    class StaticOcpQpTest 
    :   public Test
    {
    public:
        using Real = double;
        static size_t constexpr NX = 4;
        static size_t constexpr NU = 3;
        static size_t constexpr NC = 2;

        using Qp = StaticOcpQp<Real, NX, NU, NC>;

    protected:
        StaticOcpQpTest()
        :   graph_(3)
        {
        }


        OcpTree const graph_;
        Qp qp_ {graph_};
    };


    TEST_F(StaticOcpQpTest, testQpPropertyMapInterface)
    {
        auto& g = this->qp_.graph();
        auto const N = num_vertices(g);

        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>>> Q;
        std::map<OcpVertex, blaze::StaticVector<Real, NX>> q;
        std::map<OcpVertex, blaze::StaticVector<Real, NU>> r;
        std::map<OcpVertex, blaze::StaticMatrix<Real, NU, NX, blaze::columnMajor>> S;
        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NU, NU, blaze::columnMajor>>> R;

        std::map<OcpEdge, blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>> A;
        std::map<OcpEdge, blaze::StaticMatrix<Real, NX, NU, blaze::columnMajor>> B;
        std::map<OcpEdge, blaze::StaticVector<Real, NX>> b;

        std::map<OcpVertex, blaze::StaticVector<Real, NX>> x_min, x_max;
        std::map<OcpVertex, blaze::StaticVector<Real, NU>> u_min, u_max;
        std::map<OcpVertex, blaze::StaticVector<Real, NC>> d_min, d_max;

        // Writing random data
        for (auto v : vertices(g))
        {
            DynamicOcpSize const sz {NX, NU, NC};

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


    TEST_F(StaticOcpQpTest, testQpInterface)
    {
        auto& g = this->qp_.graph();
        auto const N = num_vertices(g);

        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>>> Q;
        std::map<OcpVertex, blaze::StaticVector<Real, NX>> q;
        std::map<OcpVertex, blaze::StaticVector<Real, NU>> r;
        std::map<OcpVertex, blaze::StaticMatrix<Real, NU, NX, blaze::columnMajor>> S;
        std::map<OcpVertex, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NU, NU, blaze::columnMajor>>> R;

        std::map<OcpEdge, blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>> A;
        std::map<OcpEdge, blaze::StaticMatrix<Real, NX, NU, blaze::columnMajor>> B;
        std::map<OcpEdge, blaze::StaticVector<Real, NX>> b;

        std::map<OcpVertex, blaze::StaticVector<Real, NX>> x_min, x_max;
        std::map<OcpVertex, blaze::StaticVector<Real, NU>> u_min, u_max;
        std::map<OcpVertex, blaze::StaticVector<Real, NC>> d_min, d_max;

        // Writing random data
        for (auto v : vertices(g))
        {
            DynamicOcpSize const sz {NX, NU, NC};

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


    TEST_F(StaticOcpQpTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->qp_;

        for (auto v : vertices(ws.graph()))
        {
            EXPECT_EQ(rows(this->qp_.Q(v)), NX);
            EXPECT_EQ(columns(this->qp_.Q(v)), NX);
            EXPECT_EQ(size(this->qp_.q(v)), NX);            
            EXPECT_EQ(size(this->qp_.lx(v)), NX);
            EXPECT_EQ(size(this->qp_.ux(v)), NX);            
            EXPECT_EQ(size(this->qp_.ld(v)), NC);
            EXPECT_EQ(size(this->qp_.ud(v)), NC);
        }

        for (auto v : ws.graph().branchVertices())
        {
            EXPECT_EQ(rows(this->qp_.R(v)), NU);
            EXPECT_EQ(columns(this->qp_.R(v)), NU);
            EXPECT_EQ(rows(this->qp_.S(v)), NU);
            EXPECT_EQ(columns(this->qp_.S(v)), NX);
            EXPECT_EQ(size(this->qp_.r(v)), NU);
            EXPECT_EQ(size(this->qp_.lu(v)), NU);
            EXPECT_EQ(size(this->qp_.uu(v)), NU);
        }

        for (auto e : edges(ws.graph()))
        {
            EXPECT_EQ(rows(this->qp_.A(e)), NX);
            EXPECT_EQ(columns(this->qp_.A(e)), NX);
            EXPECT_EQ(rows(this->qp_.B(e)), NX);
            EXPECT_EQ(columns(this->qp_.B(e)), NU);
            EXPECT_EQ(size(this->qp_.b(e)), NX);
        }
    }
}
