/*
* Unit tests for common functionality of StaticOcpQp class.
*/

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
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
        {
        }


        OcpGraph const graph_ = ocpGraphLinear(3);
        Qp qp_ {graph_};
    };


    TEST_F(StaticOcpQpTest, testQpInterface)
    {
        auto& g = this->qp_.graph();
        auto const N = num_vertices(g);

        std::map<OcpVertexDescriptor, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>>> Q;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NX>> q;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NU>> r;
        std::map<OcpVertexDescriptor, blaze::StaticMatrix<Real, NU, NX, blaze::columnMajor>> S;
        std::map<OcpVertexDescriptor, blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NU, NU, blaze::columnMajor>>> R;

        std::map<OcpEdgeDescriptor, blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>> A;
        std::map<OcpEdgeDescriptor, blaze::StaticMatrix<Real, NX, NU, blaze::columnMajor>> B;
        std::map<OcpEdgeDescriptor, blaze::StaticVector<Real, NX>> b;

        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NX>> x_min, x_max;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NU>> u_min, u_max;
        std::map<OcpVertexDescriptor, blaze::StaticVector<Real, NC>> d_min, d_max;

        // Writing random data
        for (auto v : graph::vertices(g))
        {
            OcpSize const sz {NX, NU, NC};

            makePositiveDefinite(Q[v]);
            makePositiveDefinite(R[v]);
            randomize(S[v]);
            randomize(q[v]);
            randomize(r[v]);
            randomize(x_min[v]);
            randomize(x_max[v]);
            randomize(u_min[v]);
            randomize(u_max[v]);
            randomize(d_min[v]);
            randomize(d_max[v]);

            put(this->qp_.Q(), v, Q[v]);
            put(this->qp_.R(), v, R[v]);
            put(this->qp_.S(), v, S[v]);
            put(this->qp_.q(), v, q[v]);
            put(this->qp_.r(), v, r[v]);

            put(this->qp_.lx(), v, x_min[v]);
            put(this->qp_.ux(), v, x_max[v]);
            put(this->qp_.lu(), v, u_min[v]);
            put(this->qp_.uu(), v, u_max[v]);
            put(this->qp_.ld(), v, d_min[v]);
            put(this->qp_.ud(), v, d_max[v]);
        }

        for (auto e : graph::edges(g))
        {
            randomize(A[e]);
            randomize(B[e]);
            randomize(b[e]);

            put(this->qp_.A(), e, A[e]);
            put(this->qp_.B(), e, B[e]);
            put(this->qp_.b(), e, b[e]);
        }

        // Reading the data and checking that they are the same that we wrote
        for (auto v : graph::vertices(g))
        {
            EXPECT_EQ(get(qp_.size(), v), OcpSize(NX, NU, NC));

            EXPECT_EQ(forcePrint(get(this->qp_.Q(), v)), forcePrint(Q[v])) << "at v=" << v;
            
            EXPECT_EQ(forcePrint(get(this->qp_.R(), v)), forcePrint(R[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.S(), v)), forcePrint(S[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.q(), v)), forcePrint(q[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.r(), v)), forcePrint(r[v])) << "at v=" << v;

            EXPECT_EQ(forcePrint(get(this->qp_.lx(), v)), forcePrint(x_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.ux(), v)), forcePrint(x_max[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.lu(), v)), forcePrint(u_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.uu(), v)), forcePrint(u_max[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.ld(), v)), forcePrint(d_min[v])) << "at v=" << v;
            EXPECT_EQ(forcePrint(get(this->qp_.ud(), v)), forcePrint(d_max[v])) << "at v=" << v;
        }

        for (auto e : graph::edges(g))
        {
            EXPECT_EQ(forcePrint(get(this->qp_.A(), e)), forcePrint(A[e])) << "at e=" << get(graph::edge_index, g, e);
            EXPECT_EQ(forcePrint(get(this->qp_.B(), e)), forcePrint(B[e])) << "at e=" << get(graph::edge_index, g, e);
            EXPECT_EQ(forcePrint(get(this->qp_.b(), e)), forcePrint(b[e])) << "at e=" << get(graph::edge_index, g, e);
        }
    }


    TEST_F(StaticOcpQpTest, testMatrixSizesCorrect)
    {
        auto const& ws = this->qp_;

        for (auto v : graph::vertices(ws.graph()))
        {
            OcpSize const s {NX, NU, NC};
            
            EXPECT_EQ(rows(get(this->qp_.Q(), v)), s.nx());
            EXPECT_EQ(columns(get(this->qp_.Q(), v)), s.nx());

            EXPECT_EQ(rows(get(this->qp_.R(), v)), s.nu());
            EXPECT_EQ(columns(get(this->qp_.R(), v)), s.nu());
            EXPECT_EQ(rows(get(this->qp_.S(), v)), s.nu());
            EXPECT_EQ(columns(get(this->qp_.S(), v)), s.nx());
            EXPECT_EQ(size(get(this->qp_.q(), v)), s.nx());
            EXPECT_EQ(size(get(this->qp_.r(), v)), s.nu());

            EXPECT_EQ(size(get(this->qp_.lx(), v)), s.nx());
            EXPECT_EQ(size(get(this->qp_.ux(), v)), s.nx());
            EXPECT_EQ(size(get(this->qp_.lu(), v)), s.nu());
            EXPECT_EQ(size(get(this->qp_.uu(), v)), s.nu());
            EXPECT_EQ(size(get(this->qp_.ld(), v)), s.nc());
            EXPECT_EQ(size(get(this->qp_.ud(), v)), s.nc());
        }

        for (auto e : graph::edges(ws.graph()))
        {
            OcpSize const s {NX, NU, NC};

            EXPECT_EQ(rows(get(this->qp_.A(), e)), s.nx());
            EXPECT_EQ(columns(get(this->qp_.A(), e)), s.nx());
            EXPECT_EQ(rows(get(this->qp_.B(), e)), s.nx());
            EXPECT_EQ(columns(get(this->qp_.B(), e)), s.nu());
            EXPECT_EQ(size(get(this->qp_.b(), e)), s.nx());
        }
    }
}
