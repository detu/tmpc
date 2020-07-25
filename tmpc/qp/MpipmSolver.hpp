#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>

#include <vector>


namespace tmpc
{
    template <typename Real_>
    class MpipmSolver
    {
    public:
        using Real = Real_;

        
        MpipmSolver(DynamicOcpSize const& size)
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


        auto const& size() const noexcept
        {
            return size_;
        }


        auto const& graph() const
        {
            return graph_;
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
            ,   x_(sz.nx(v))
            ,   u_(sz.nu(v))
            ,   lam_lx_(sz.nx(v))
            ,   lam_ux_(sz.nx(v))
            ,   lam_lu_(sz.nu(v))
            ,   lam_uu_(sz.nu(v))
            ,   lam_ld_(sz.nc(v))
            ,   lam_ud_(sz.nc(v))
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

            blaze::DynamicVector<Real> x_;
            blaze::DynamicVector<Real> u_;
            blaze::DynamicVector<Real> lam_lx_;
            blaze::DynamicVector<Real> lam_ux_;
            blaze::DynamicVector<Real> lam_lu_;
            blaze::DynamicVector<Real> lam_uu_;
            blaze::DynamicVector<Real> lam_ld_;
            blaze::DynamicVector<Real> lam_ud_;
        };


        struct EdgePropertyBundle
        {
            EdgePropertyBundle(DynamicOcpSize const& size, OcpEdge e)
            :   A_(size.nx(target(e, size.graph())), size.nx(source(e, size.graph())))
            ,   B_(size.nx(target(e, size.graph())), size.nu(source(e, size.graph())))
			/*
            :   BA_(size_dst.nx(), size_src.nu() + size_src.nx())
            ,   A_(submatrix(BA_, 0, size_src.nu(), size_dst.nx(), size_src.nx()))
            ,   B_(submatrix(BA_, 0, 0, size_dst.nx(), size_src.nu()))
			*/
            ,   b_(size.nx(target(e, size.graph())))
            ,   pi_(size.nx(target(e, size.graph())))
            {
            }


            blaze::DynamicMatrix<Real, blaze::columnMajor> A_;
            blaze::DynamicMatrix<Real, blaze::columnMajor> B_;
			
			/*
            blaze::DynamicMatrix<Real, blaze::columnMajor> BA_;

            using Submatrix = decltype(blaze::submatrix(BA_, 0, 0, 1, 1));

            Submatrix A_;
            Submatrix B_;
			*/
			
            blaze::DynamicVector<Real> b_;
            blaze::DynamicVector<Real> pi_;
        };


        OcpTree graph_;
        DynamicOcpSize size_;
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };
}