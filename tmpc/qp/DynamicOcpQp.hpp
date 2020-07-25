#pragma once

#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>

#include <vector>
#include <initializer_list>


namespace tmpc
{
    template <typename Real_>
    class DynamicOcpQp
    {
    public:
        using Real = Real_;


        DynamicOcpQp()
        :   DynamicOcpQp(DynamicOcpSize {{0, 0, 0}})
        {
        }

        
        template <OcpSize Size>
        explicit DynamicOcpQp(Size const& size)
        :   graph_ {size.graph()}
        ,   size_ {size}
        {
            // Allocate vertex properties of appropriate size
            vertexProperties_.reserve(num_vertices(graph_));
            for (auto v : vertices(graph_))
                vertexProperties_.emplace_back(size_, v);

            // Populate edgeIndex_ and allocate edge properties of appropriate size
            edgeProperties_.reserve(num_edges(graph_));
            for (auto e : edges(graph_))
                edgeProperties_.emplace_back(size_, e);
        }


        /// @brief OCP size.
        auto const& size() const
        {
            return size_;
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto const& H(OcpVertex v) const
        {
            return vertexProperties_[v].H_;
        }


        template <typename MT, bool SO>
        void H(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            noresize(vertexProperties_[v].H_) = ~val;
        }


        auto const& Q(OcpVertex v) const
        {
            return vertexProperties_[v].Q_;
        }


        auto& Q(OcpVertex v)
        {
            return vertexProperties_[v].Q_;
        }


        template <typename MT, bool SO>
        void Q(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            vertexProperties_[v].Q_ = ~val;
        }


        auto const& R(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].R_;
        }


        template <typename MT, bool SO>
        void R(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            vertexProperties_[v].R_ = val;
        }


        template <typename MT, bool SO>
        void S(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            vertexProperties_[v].S_ = ~val;
        }


        auto const& S(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].S_;
        }


        template <typename VT, bool TF>
        void q(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            noresize(vertexProperties_[v].q_) = ~val;
        }


        auto const& q(OcpVertex v) const
        {
            return vertexProperties_[v].q_;
        }


        template <typename VT, bool TF>
        void r(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            noresize(vertexProperties_[v].r_) = ~val;
        }


        auto const& r(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].r_;
        }


        template <typename MT, bool SO>
        void C(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            noresize(vertexProperties_[v].C_) = ~val;
        }


        auto const& C(OcpVertex v) const
        {
            return vertexProperties_[v].C_;
        }


        template <typename MT, bool SO>
        void D(OcpVertex v, blaze::Matrix<MT, SO> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            noresize(vertexProperties_[v].D_) = ~val;
        }


        auto const& D(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].D_;
        }


        template <typename VT, bool TF>
        void ld(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            noresize(vertexProperties_[v].ld_) = ~val;
        }


        auto const& ld(OcpVertex v) const
        {
            return vertexProperties_[v].ld_;
        }


        template <typename VT, bool TF>
        void ud(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            noresize(vertexProperties_[v].ud_) = ~val;
        }


        auto const& ud(OcpVertex v) const
        {
            return vertexProperties_[v].ud_;
        }


        template <typename VT, bool TF>
        void lx(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            noresize(vertexProperties_[v].lx_) = ~val;
        }


        auto const& lx(OcpVertex v) const
        {
            return vertexProperties_[v].lx_;
        }


        template <typename VT, bool TF>
        void ux(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            noresize(vertexProperties_[v].ux_) = ~val;
        }


        auto const& ux(OcpVertex v) const
        {
            return vertexProperties_[v].ux_;
        }


        template <typename VT, bool TF>
        void lu(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            noresize(vertexProperties_[v].lu_) = ~val;
        }


        auto const& lu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].lu_;
        }


        template <typename VT, bool TF>
        void uu(OcpVertex v, blaze::Vector<VT, TF> const& val)
        {
            assert(out_degree(v, graph_) > 0);
            noresize(vertexProperties_[v].uu_) = ~val;
        }


        auto const& uu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].uu_;
        }


        auto const& BA(OcpEdge e) const
        {
            return edgeProperties_[e].BA_;
        }


        template <typename MT, bool SO>
        void BA(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            noresize(edgeProperties_[e].BA_) = ~val;
        }


        auto const& A(OcpEdge e) const
        {
            return edgeProperties_[e].A_;
        }


        template <typename MT, bool SO>
        void A(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            edgeProperties_[e].A_ = ~val;
        }


        auto const& B(OcpEdge e) const
        {
            return edgeProperties_[e].B_;
        }


        template <typename MT, bool SO>
        void B(OcpEdge e, blaze::Matrix<MT, SO> const& val)
        {
            edgeProperties_[e].B_ = ~val;
        }


        auto const& b(OcpEdge e) const
        {
            return edgeProperties_[e].b_;
        }


        template <typename VT, bool TF>
        void b(OcpEdge e, blaze::Vector<VT, TF> const& val)
        {
            noresize(edgeProperties_[e].b_) = ~val;
        }


    private:
        struct VertexPropertyBundle
        {
            VertexPropertyBundle(DynamicOcpSize const& sz, OcpVertex v)
            :   H_(sz.nu(v) + sz.nx(v))
            ,   Q_(submatrix(H_, sz.nu(v), sz.nu(v), sz.nx(v), sz.nx(v)))
            ,   R_(submatrix(H_, 0, 0, sz.nu(v), sz.nu(v)))
            ,   S_(submatrix(H_, 0, sz.nu(v), sz.nu(v), sz.nx(v)))
            ,   q_(sz.nx(v))
            ,   r_(sz.nu(v))
            ,   lx_(sz.nx(v))
            ,   ux_(sz.nx(v))
            ,   lu_(sz.nu(v))
            ,   uu_(sz.nu(v))
            ,   C_(sz.nc(v), sz.nx(v))
            ,   D_(sz.nc(v), sz.nu(v))
            ,   ld_(sz.nc(v))
            ,   ud_(sz.nc(v))
            {
            }

            // H = [ R,   S
            //       S^T, Q]
            blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> H_;

            using Submatrix = decltype(blaze::submatrix(H_, 0, 0, 1, 1));

            Submatrix Q_;
            Submatrix R_;
            Submatrix S_;
            blaze::DynamicVector<Real> q_;
            blaze::DynamicVector<Real> r_;

            blaze::DynamicVector<Real> lx_;
            blaze::DynamicVector<Real> ux_;
            blaze::DynamicVector<Real> lu_;
            blaze::DynamicVector<Real> uu_;

            blaze::DynamicMatrix<Real, blaze::columnMajor> C_;
            blaze::DynamicMatrix<Real, blaze::columnMajor> D_;
            blaze::DynamicVector<Real> ld_;
            blaze::DynamicVector<Real> ud_;
        };


        struct EdgePropertyBundle
        {
            EdgePropertyBundle(DynamicOcpSize const& size, OcpEdge e)
            :   BA_(targetNx(size, e), sourceNu(size, e) + sourceNx(size, e))
			,   A_(submatrix(BA_, 0, sourceNu(size, e), targetNx(size, e), sourceNx(size, e)))
            ,   B_(submatrix(BA_, 0, 0, targetNx(size, e), sourceNu(size, e)))
			,   b_(targetNx(size, e))
            {
            }


            blaze::DynamicMatrix<Real, blaze::columnMajor> BA_;

            using Submatrix = decltype(blaze::submatrix(BA_, 0, 0, 1, 1));

            Submatrix A_;
            Submatrix B_;
			
            blaze::DynamicVector<Real> b_;
        };


        OcpTree graph_;
        DynamicOcpSize size_;
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };
}