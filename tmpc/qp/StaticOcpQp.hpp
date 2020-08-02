#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    template <typename Real_, size_t NX, size_t NU, size_t NC = 0>
    class StaticOcpQp
    {
    public:
        using Real = Real_;


        StaticOcpQp()
        :   graph_ {}
        ,   size_ {graph_}
        ,   edgeVertexProperties_(num_vertices(graph_))
        {
        }


        explicit StaticOcpQp(OcpTree const& g)
        :   graph_(g)
        ,   size_ {graph_}
        ,   edgeVertexProperties_(num_vertices(g))
        {
        }


        template <OcpQp Qp>
        StaticOcpQp& operator=(Qp const& rhs)
        {
            assign(*this, rhs);
            return *this;
        }


        /// @brief OCP size.
        auto const& size() const noexcept
        {
            return size_;
        }


        auto const& graph() const noexcept
        {
            return graph_;
        }


        auto const& H(OcpVertex v) const
        {
            return edgeVertexProperties_[v].H_;
        }


        template <typename MT, bool SO>
        void H(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            edgeVertexProperties_[v].H_ = ~val;
        }


        auto const& Q(OcpVertex v) const
        {
            return edgeVertexProperties_[v].Q_;
        }


        template <typename MT, bool SO>
        void Q(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            edgeVertexProperties_[v].Q_ = ~val;
        }


        auto const& R(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].R_;
        }


        template <typename MT, bool SO>
        void R(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].R_ = val;
        }


        template <typename MT, bool SO>
        void S(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].S_ = ~val;
        }


        auto const& S(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].S_;
        }


        template <typename VT, bool TF>
        void q(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[v].q_ = ~val;
        }


        auto const& q(OcpVertex v) const
        {
            return edgeVertexProperties_[v].q_;
        }


        template <typename VT, bool TF>
        void r(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].r_ = ~val;
        }


        auto const& r(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].r_;
        }


        template <typename MT, bool SO>
        void C(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            edgeVertexProperties_[v].C_ = ~val;
        }


        auto const& C(OcpVertex v) const
        {
            return edgeVertexProperties_[v].C_;
        }


        template <typename MT, bool SO>
        void D(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].D_ = ~val;
        }


        auto const& D(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].D_;
        }


        template <typename VT, bool TF>
        void ld(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[v].ld_ = ~val;
        }


        auto const& ld(OcpVertex v) const
        {
            return edgeVertexProperties_[v].ld_;
        }


        template <typename VT, bool TF>
        void ud(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[v].ud_ = ~val;
        }


        auto const& ud(OcpVertex v) const
        {
            return edgeVertexProperties_[v].ud_;
        }


        template <typename VT, bool TF>
        void lx(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[v].lx_ = ~val;
        }


        auto const& lx(OcpVertex v) const
        {
            return edgeVertexProperties_[v].lx_;
        }


        template <typename VT, bool TF>
        void ux(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[v].ux_ = ~val;
        }


        auto const& ux(OcpVertex v) const
        {
            return edgeVertexProperties_[v].ux_;
        }


        template <typename VT, bool TF>
        void lu(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].lu_ = ~val;
        }


        auto const& lu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].lu_;
        }


        template <typename VT, bool TF>
        void uu(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            edgeVertexProperties_[v].uu_ = ~val;
        }


        auto const& uu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return edgeVertexProperties_[v].uu_;
        }


        auto const& BA(OcpEdge e) const
        {
            return edgeVertexProperties_[target(e, graph_)].BA_;
        }


        template <typename MT, bool SO>
        void BA(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            edgeVertexProperties_[target(e, graph_)].BA_ = ~val;
        }


        auto A(OcpEdge e) const
        {
            return blaze::submatrix<0, NU, NX, NX>(edgeVertexProperties_[target(e, graph_)].BA_);
        }


        template <typename MT, bool SO>
        void A(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            blaze::submatrix<0, NU, NX, NX>(edgeVertexProperties_[target(e, graph_)].BA_) = ~val;
        }


        auto B(OcpEdge e) const
        {
            return blaze::submatrix<0, 0, NX, NU>(edgeVertexProperties_[target(e, graph_)].BA_);
        }


        template <typename MT, bool SO>
        void B(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            blaze::submatrix<0, 0, NX, NU>(edgeVertexProperties_[target(e, graph_)].BA_) = ~val;
        }


        auto const& b(OcpEdge e) const
        {
            return edgeVertexProperties_[target(e, graph_)].b_;
        }


        template <typename VT, bool TF>
        void b(OcpEdge e, blaze::Vector<VT, TF> const& val)
        {
            edgeVertexProperties_[target(e, graph_)].b_ = ~val;
        }


    private:
        /// @brief Stores properties of a vertex and its parent edge together to improve data locality.
        struct EdgeVertexPropertyBundle
        {
            //
            // Edge properties
            //

            blaze::StaticMatrix<Real, NX, NU + NX, blaze::rowMajor> BA_;
			blaze::StaticVector<Real, NX> b_;

            //
            // Vertex properties
            //

            // H = [ R,   S
            //       S^T, Q]
            blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NU + NX, NU + NX, blaze::columnMajor>> H_;

            decltype(blaze::submatrix<NU, NU, NX, NX>(H_)) Q_ = blaze::submatrix<NU, NU, NX, NX>(H_);
            decltype(blaze::submatrix<0, 0, NU, NU>(H_)) R_ = blaze::submatrix<0, 0, NU, NU>(H_);
            decltype(blaze::submatrix<0, NU, NU, NX>(H_)) S_ = blaze::submatrix<0, NU, NU, NX>(H_);
            blaze::StaticVector<Real, NX> q_;
            blaze::StaticVector<Real, NU> r_;

            blaze::StaticVector<Real, NX> lx_;
            blaze::StaticVector<Real, NX> ux_;
            blaze::StaticVector<Real, NU> lu_;
            blaze::StaticVector<Real, NU> uu_;

            blaze::StaticMatrix<Real, NC, NX, blaze::columnMajor> C_;
            blaze::StaticMatrix<Real, NC, NU, blaze::columnMajor> D_;
            blaze::StaticVector<Real, NC> ld_;
            blaze::StaticVector<Real, NC> ud_;
        };


        OcpTree graph_;
        StaticOcpSize<NX, NU, NC> size_;
        std::vector<EdgeVertexPropertyBundle> edgeVertexProperties_;
    };
}