#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/qp/detail/MatrixPropertyMap.hpp>
#include <tmpc/qp/detail/VectorPropertyMap.hpp>

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
            opts_.lineSearchBeta = 0.8;
            return *this;
        }


        treeqp_tdunes_opts_t const& nativeOptions() const
        {
            return opts_;
        }


    private:
        treeqp_tdunes_opts_t opts_;
        std::unique_ptr<char []> mem_;
    };


    /// Visitor that writes out-degree of each vertex to an output iterator.
    /// It also checks for non-tree edges and throws an std::invalid_argument() if one is found.
    ///
    template <typename OutIter>
    class OutDegreeVisitor
    :   public graph::default_bfs_visitor 
    {
    public:
        OutDegreeVisitor(OutIter iter)
        :   iter_{iter}
        {
        }

    
        template <typename Vertex, typename Graph>
        void discover_vertex(Vertex u, const Graph& g)
        {
            *iter_++ = out_degree(u, g);
        }

    
    private:
        OutIter iter_;
    };


    template <typename Kernel_>
    class DualNewtonTreeWorkspace
    {
    public:
        using Kernel = Kernel_;
        using Real = typename Kernel::Real;


        template <typename SizeMap>
        DualNewtonTreeWorkspace(OcpGraph const& g, SizeMap sz)
        :   DualNewtonTreeWorkspace(g, sz, DualNewtonTreeOptions(num_vertices(g)))
        {
        }


        template <typename SizeMap>
        DualNewtonTreeWorkspace(OcpGraph const& g, SizeMap sz, DualNewtonTreeOptions&& options)
        :   graph_{g}
        ,   size_{num_vertices(g)}
        ,   tree_{num_vertices(g)}
        ,   opts_(std::move(options))
        {
            auto const num_nodes = num_vertices(g);
            auto const vertex_id = get(graph::vertex_index, g);

            // Fill the own size_ array
            copyProperty(sz, iterator_property_map(size_.begin(), vertex_id), graph::vertices(g));

            // Fill size arrays.
            std::vector<int> nx(num_nodes), nu(num_nodes), nc(num_nodes), nk(num_nodes);

            for (auto v : graph::vertices(g))
            {
                auto const i = vertex_id[v];
                nx[i] = size_[i].nx();
                nu[i] = size_[i].nu();
                nc[i] = size_[i].nc();
            }
            
            // Fill the number of kids vector with out-degrees of the nodes.
            breadth_first_search(g, vertex(0, g), visitor(OutDegreeVisitor {nk.begin()}));

            // Setup the tree.
            auto const qp_in_size = tree_qp_in_calculate_size(
                num_nodes, nx.data(), nu.data(), nc.data(), nk.data());

            qp_in_memory_.resize(qp_in_size);
            tree_qp_in_create(num_nodes, nx.data(), nu.data(), nc.data(), nk.data(), 
                &qp_in_, qp_in_memory_.data());

            auto const qp_out_size = tree_qp_out_calculate_size(
                num_nodes, nx.data(), nu.data(), nc.data());

            qp_out_memory_.resize(qp_out_size);
            tree_qp_out_create(num_nodes, nx.data(), nu.data(), nc.data(),
                &qp_out_, qp_out_memory_.data());

            auto const treeqp_size = treeqp_tdunes_calculate_size(&qp_in_, &opts_.nativeOptions());
            qp_solver_memory_.resize(treeqp_size);
            treeqp_tdunes_create(&qp_in_, &opts_.nativeOptions(), &work_, qp_solver_memory_.data());
        }


        void solve()
        {
            // A workaround for this issue: https://gitlab.syscop.de/dimitris.kouzoupis/hangover/issues/6
            if (true)
            {
                size_t sum_nx = 0;
                for (auto v : graph::vertices(graph_))
                    sum_nx += get(size(), v).nx();

                std::vector<Real> lambda(sum_nx);
                treeqp_tdunes_set_dual_initialization(lambda.data(), &work_);
            }
            
            auto const ret = treeqp_tdunes_solve(&qp_in_, &qp_out_, &opts_.nativeOptions(), &work_);

            if (ret != TREEQP_OPTIMAL_SOLUTION_FOUND)
                throw DualNewtonTreeException(ret);
        }


        void print() const
        {
            tree_qp_in_print_dims(&qp_in_);
            tree_qp_in_print(&qp_in_);
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto edgeIndex() const
        {
            return make_function_property_map<OcpEdgeDescriptor>(
                [this] (OcpEdgeDescriptor e)
                {
                    return target(e, graph_) - 1;
                }
            );
        }


        auto size() const
        {
            return iterator_property_map(size_.begin(), vertexIndex(graph_));
        }


        auto Q()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_Q(size()), 
                std::bind(tree_qp_in_set_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto Q() const
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_Q(size()), 
                std::bind(tree_qp_in_set_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto R()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_R(size()), 
                std::bind(tree_qp_in_set_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto R() const
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_R(size()), 
                std::bind(tree_qp_in_set_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto S()
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc now assumes the same!
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_S(size()), 
                std::bind(tree_qp_in_set_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto S() const
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc now assumes the same!
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_S(size()), 
                std::bind(tree_qp_in_set_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto q()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto q() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto r()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto r() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto C()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_C(size()), 
                std::bind(tree_qp_in_set_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto C() const
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_C(size()), 
                std::bind(tree_qp_in_set_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto D()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_D(size()), 
                std::bind(tree_qp_in_set_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto D() const
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(graph_), size_D(size()), 
                std::bind(tree_qp_in_set_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto ld()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_d(size()), 
                std::bind(tree_qp_in_set_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ld() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_d(size()), 
                std::bind(tree_qp_in_set_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ud()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_d(size()), 
                std::bind(tree_qp_in_set_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ud() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_d(size()), 
                std::bind(tree_qp_in_set_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lx()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lx() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ux()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ux() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_in_set_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lu()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lu() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto uu()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto uu() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_in_set_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto A()
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_A(size(), graph_), 
                std::bind(tree_qp_in_set_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto A() const
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_A(size(), graph_), 
                std::bind(tree_qp_in_set_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto B()
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_B(size(), graph_), 
                std::bind(tree_qp_in_set_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto B() const
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_B(size(), graph_), 
                std::bind(tree_qp_in_set_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_qp_in_get_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto b()
        {
            // void tree_qp_in_set_edge_b(const double * const b, tree_qp_in * const qp_in, const int indx);
            // void tree_qp_in_get_edge_b(double * const b, const tree_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpEdgeDescriptor, DynamicVector<Kernel>>(edgeIndex(), size_b(size(), graph_), 
                std::bind(tree_qp_in_set_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto b() const
        {
            // void tree_qp_in_set_edge_b(const double * const b, tree_qp_in * const qp_in, const int indx);
            // void tree_qp_in_get_edge_b(double * const b, const tree_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpEdgeDescriptor, DynamicVector<Kernel>>(edgeIndex(), size_b(size(), graph_), 
                std::bind(tree_qp_in_set_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_qp_in_get_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto x() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_out_set_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto u() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_out_set_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_x() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_x(size()), 
                std::bind(tree_qp_out_set_node_mu_x, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_x, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_u() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_u(size()), 
                std::bind(tree_qp_out_set_node_mu_u, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_u, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto mu_d() const
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(graph_), size_d(size()), 
                std::bind(tree_qp_out_set_node_mu_d, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_node_mu_d, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto pi() const
        {
            return detail::makeVectorPropertyMap<OcpEdgeDescriptor, DynamicVector<Kernel>>(get(graph::edge_index, graph_), size_b(size(), graph_), 
                std::bind(tree_qp_out_set_edge_lam, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_qp_out_get_edge_lam, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }
        

    private:
        OcpGraph graph_;
        std::vector<OcpSize> size_;

        std::vector<node> tree_;
        tree_qp_in qp_in_;
        std::vector<char> qp_in_memory_;
        tree_qp_out qp_out_;
        std::vector<char> qp_out_memory_;
        DualNewtonTreeOptions opts_;
        treeqp_tdunes_workspace work_;
        std::vector<char> qp_solver_memory_;
    };
}