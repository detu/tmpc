#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/property_map/MatrixPropertyMap.hpp>
#include <tmpc/property_map/VectorPropertyMap.hpp>

#include <tmpc/Traits.hpp>

#define TMPC_TREEQP_USE_HPMPC

#ifdef TMPC_TREEQP_USE_HPMPC
#include <treeqp/src/hpmpc_tree.h>
#endif

#include <treeqp/src/tree_qp_common.h>
#include <treeqp/src/dual_Newton_tree.h>
#include <treeqp/utils/tree.h>
#include <treeqp/utils/print.h>

#include <vector>
#include <memory>
#include <stdexcept>
#include <functional>


namespace tmpc
{
    class DualNewtonTreeException
	:	public std::runtime_error
	{
	public:
		DualNewtonTreeException(return_t code);
		return_t code() const { return _code; }

	private:
		return_t const _code;
	};


    /// @brief Encapsulates options for the treeQP Dual Newton solver.
    class DualNewtonTreeOptions
    {
    public:
        DualNewtonTreeOptions(size_t num_nodes);


        /// @brief Move ctor
        DualNewtonTreeOptions(DualNewtonTreeOptions&&) = default;


        DualNewtonTreeOptions& maxIter(size_t max_iter)
        {
            opts_.maxIter = max_iter;
            return *this;
        }

        #ifndef TMPC_TREEQP_USE_HPMPC

        DualNewtonTreeOptions& stationarityTolerance(double val)
        {
            opts_.stationarityTolerance = val;
            return *this;
        }


        DualNewtonTreeOptions& lineSearchMaxIter(size_t val)
        {
            opts_.lineSearchMaxIter = val;
            return *this;
        }


        DualNewtonTreeOptions& lineSearchBeta(double val)
        {
            opts_.lineSearchBeta = val;
            return *this;
        }


        DualNewtonTreeOptions& lineSearchGamma(double val)
        {
            opts_.lineSearchGamma = val;
            return *this;
        }


        DualNewtonTreeOptions& lineSearchTol(double val)
        {
            opts_.lineSearchTol = val;
            return *this;
        }
		
		
        treeqp_tdunes_opts_t const& nativeOptions() const
        {
            return opts_;
        }


        DualNewtonTreeOptions& lineSearchRestartTrigger(int val)
        {
            opts_.lineSearchRestartTrigger = val;
            return *this;
        }


        DualNewtonTreeOptions& regType(std::string val)
        {
            if (val == "TREEQP_ALWAYS_LEVENBERG_MARQUARDT")
                opts_.regType = TREEQP_ALWAYS_LEVENBERG_MARQUARDT;
            else if (val == "TREEQP_ON_THE_FLY_LEVENBERG_MARQUARDT")
                opts_.regType = TREEQP_ON_THE_FLY_LEVENBERG_MARQUARDT;
            else
                opts_.regType = TREEQP_NO_REGULARIZATION;

            return *this;
        }


        DualNewtonTreeOptions& regTol(double val)
        {
            opts_.regTol = val;
            return *this;
        }


        DualNewtonTreeOptions& regValue(double val)
        {
            opts_.regValue = val;
            return *this;
        }

        #endif

        #ifdef TMPC_TREEQP_USE_HPMPC
        treeqp_hpmpc_opts_t const& nativeOptions() const
        {
            return opts_;
        }
        #else
        treeqp_tdunes_opts_t const& nativeOptions() const
        {
            return opts_;
        }
        #endif

    private:

        #ifdef TMPC_TREEQP_USE_HPMPC
        treeqp_hpmpc_opts_t opts_;
        #else
        treeqp_tdunes_opts_t opts_;
        #endif
        std::unique_ptr<char []> mem_;
    };


    class DualNewtonTreeWorkspace
    {
        auto edgeIndex() const
        {
            return edgeIndex(graph_);
        }


    public:
        using Real = double;


        template <typename SizeMap>
        DualNewtonTreeWorkspace(OcpTree const& g, SizeMap sz)
        :   DualNewtonTreeWorkspace(g, sz, DualNewtonTreeOptions(num_vertices(g)))
        {
        }


        template <typename SizeMap>
        DualNewtonTreeWorkspace(OcpTree const& g, SizeMap sz, DualNewtonTreeOptions&& options)
        :   graph_{g}
        ,   size_{num_vertices(g)}
        ,   tree_{num_vertices(g)}
        ,   opts_(std::move(options))
        {
            // Fill the own size_ array
            copyProperty(sz, iterator_property_map(size_.begin(), vertexIndex(g)), vertices(g));

            // Finalize construction
            init();
        }


        auto size() const
        {
            return iterator_property_map(size_.begin(), vertexIndex(graph_));
        }


        /// @brief Solve the QP
        void solve();


        /// @brief Print QP dimensions and data using treeQP functions.
        void print() const;


        auto const& graph() const
        {
            return graph_;
        }


        auto Q()
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_Q(size()),
                std::bind(tree_qp_in_set_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto Q() const
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_Q(size()),
                std::bind(tree_qp_in_set_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto R()
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_R(size()),
                std::bind(tree_qp_in_set_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto R() const
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_R(size()),
                std::bind(tree_qp_in_set_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto S()
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc now assumes the same!
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_S(size()),
                std::bind(tree_qp_in_set_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto S() const
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc now assumes the same!
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_S(size()),
                std::bind(tree_qp_in_set_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto q()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto q() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto r()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto r() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto C()
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_C(size()),
                std::bind(tree_qp_in_set_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto C() const
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_C(size()),
                std::bind(tree_qp_in_set_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto D()
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_D(size()),
                std::bind(tree_qp_in_set_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto D() const
        {
            return detail::makeMatrixPropertyMap<OcpVertex, blaze::DynamicMatrix<Real, columnMajor>>(vertexIndex(graph_), size_D(size()),
                std::bind(tree_qp_in_set_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto ld()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_d(size()),
                std::bind(tree_qp_in_set_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ld() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_d(size()),
                std::bind(tree_qp_in_set_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ud()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_d(size()),
                std::bind(tree_qp_in_set_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ud() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_d(size()),
                std::bind(tree_qp_in_set_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lx()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lx() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ux()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ux() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_in_set_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lu()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lu() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto uu()
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto uu() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_in_set_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto A()
        {
            return detail::makeMatrixPropertyMap<OcpEdge, blaze::DynamicMatrix<Real, columnMajor>>(edgeIndex(), size_A(size(), graph_),
                std::bind(tree_qp_in_set_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto A() const
        {
            return detail::makeMatrixPropertyMap<OcpEdge, blaze::DynamicMatrix<Real, columnMajor>>(edgeIndex(), size_A(size(), graph_),
                std::bind(tree_qp_in_set_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto B()
        {
            return detail::makeMatrixPropertyMap<OcpEdge, blaze::DynamicMatrix<Real, columnMajor>>(edgeIndex(), size_B(size(), graph_),
                std::bind(tree_qp_in_set_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto B() const
        {
            return detail::makeMatrixPropertyMap<OcpEdge, blaze::DynamicMatrix<Real, columnMajor>>(edgeIndex(), size_B(size(), graph_),
                std::bind(tree_qp_in_set_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto b()
        {
            // void tree_qp_in_set_edge_b(const double * const b, tree_qp_in * const qp_in, const int indx);
            // void tree_qp_in_get_edge_b(double * const b, const tree_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpEdge, blaze::DynamicVector<Real>>(edgeIndex(), size_b(size(), graph_),
                std::bind(tree_qp_in_set_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto b() const
        {
            // void tree_qp_in_set_edge_b(const double * const b, tree_qp_in * const qp_in, const int indx);
            // void tree_qp_in_get_edge_b(double * const b, const tree_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpEdge, blaze::DynamicVector<Real>>(edgeIndex(), size_b(size(), graph_),
                std::bind(tree_qp_in_set_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto x() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_out_set_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto u() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_out_set_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_x() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_x(size()),
                std::bind(tree_qp_out_set_node_mu_x, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_x, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_u() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_u(size()),
                std::bind(tree_qp_out_set_node_mu_u, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_u, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_d() const
        {
            return detail::makeVectorPropertyMap<OcpVertex, blaze::DynamicVector<Real>>(vertexIndex(graph_), size_d(size()),
                std::bind(tree_qp_out_set_node_mu_d, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_d, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto pi() const
        {
            return detail::makeVectorPropertyMap<OcpEdge, blaze::DynamicVector<Real>>(edgeIndex(graph_), size_b(size(), graph_),
                std::bind(tree_qp_out_set_edge_lam, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_edge_lam, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


    private:
        OcpTree graph_;
        std::vector<DynamicOcpSize> size_;

        std::vector<node> tree_;
        tree_qp_in qp_in_;
        std::vector<char> qp_in_memory_;
        tree_qp_out qp_out_;
        std::vector<char> qp_out_memory_;
        DualNewtonTreeOptions opts_;
        #ifdef TMPC_TREEQP_USE_HPMPC
        treeqp_hpmpc_workspace work_;
        #else
        treeqp_tdunes_workspace work_;
        #endif
        std::vector<char> qp_solver_memory_;

        void init();
    };


    template <>
    struct RealOf<DualNewtonTreeWorkspace>
    {
        using type = double;
    };
}