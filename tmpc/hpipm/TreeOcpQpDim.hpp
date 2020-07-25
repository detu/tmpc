#pragma once

#include "Tree.hpp"

#include <hpipm_d_tree_ocp_qp_dim.h>

#include <tmpc/ocp/OcpSize.hpp>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct TreeOcpQpDimType;


    template <>
    struct TreeOcpQpDimType<double>
    {
        using Type = ::d_tree_ocp_qp_dim;
    };


    template <typename Real>
    using tree_ocp_qp_dim = typename TreeOcpQpDimType<Real>::Type;
    
    
    template <typename Real>
    inline int memsize_tree_ocp_qp_dim(int Nn);

    
    template <>
    inline int memsize_tree_ocp_qp_dim<double>(int Nn)
    {
        return ::d_memsize_tree_ocp_qp_dim(Nn);
    }


    inline void create_tree_ocp_qp_dim(int Nn, d_tree_ocp_qp_dim& qp_dim, void * memory)
    {
        ::d_create_tree_ocp_qp_dim(Nn, &qp_dim, memory);
    }


    inline void cvt_int_to_tree_ocp_qp_dim(tree const& ttree, int * nx, int * nu, 
        int * nbx, int * nbu, int * ng, int * nsbx, int * nsbu, int * nsg, d_tree_ocp_qp_dim& dim)
    {
        ::d_cvt_int_to_tree_ocp_qp_dim(const_cast<tree *>(&ttree),
            nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);
    }


    template <typename Real>
    struct TreeOcpQpDim
    :	tree_ocp_qp_dim<Real>
    {
        /// @brief Construct tree_ocp_qp_dim from OcpSize
        ///
        template <OcpSize Size>
        explicit TreeOcpQpDim(Size const& size)
        :	tree_ {size.graph()}
        ,   memory_ {new char[memsize_tree_ocp_qp_dim<Real>(num_vertices(size.graph()))]}
        {
            auto const N = num_vertices(size.graph());
            create_tree_ocp_qp_dim(N, *this, memory_.get());

            std::unique_ptr<int []> nx {new int[N]};
            std::unique_ptr<int []> nu {new int[N]};
            std::unique_ptr<int []> nbx {new int[N]};
            std::unique_ptr<int []> nbu {new int[N]};
            std::unique_ptr<int []> ng {new int[N]};
            std::unique_ptr<int []> nsbx {new int[N]};
            std::unique_ptr<int []> nsbu {new int[N]};
            std::unique_ptr<int []> nsg {new int[N]};

            for (auto v : vertices(size.graph()))
            {
                nx[v] = size.nx(v);
                nu[v] = size.nu(v);
                nbx[v] = size.nx(v);
                nbu[v] = size.nu(v);
                ng[v] = size.nc(v);
                nsbx[v] = 0;
                nsbu[v] = 0;
                nsg[v] = 0;
            }

            cvt_int_to_tree_ocp_qp_dim(tree_, nx.get(), nu.get(), nbx.get(), nbu.get(),
                ng.get(), nsbx.get(), nsbu.get(), nsg.get(), *this);
        }


    private:
        Tree tree_;
        std::unique_ptr<char []> memory_;
    };
}