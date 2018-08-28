#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/core/Graph.hpp>
#include <tmpc/qp/detail/MatrixPropertyMap.hpp>
#include <tmpc/qp/detail/VectorPropertyMap.hpp>

#include <treeqp/src/tree_ocp_qp_common.h>
#include <treeqp/src/dual_Newton_tree.h>
#include <treeqp/utils/tree.h>
#include <treeqp/utils/print.h>

#include <vector>
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


    /// Visitor that writes out-degree of each vertex to an output iterator.
    /// It also checks for non-tree edges and throws an std::invalid_argument() if one is found.
    ///
    template <typename OutIter>
    class OutDegreeVisitor
    :   public boost::default_bfs_visitor 
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


        template <typename InIterSize>
        DualNewtonTreeWorkspace(OcpGraph const& g, InIterSize sz)
        //DualNewtonTreeWorkspace(OcpGraph const& g, OcpSize const * sz)
        :   graph_{g}
        ,   size_{sz, sz + num_vertices(g)}
        ,   tree_{num_vertices(g)}
        {
            auto const num_nodes = num_vertices(g);

            // Fill size arrays.
            std::vector<int> nx(num_nodes), nu(num_nodes), nc(num_nodes), ns(num_nodes);

            for (size_t i = 0; i < num_nodes; ++i)
            {
                nx[i] = size_[i].nx();
                nu[i] = size_[i].nu();
                nc[i] = size_[i].nc();
            }
            
            // Fill the out-degree vector.
            breadth_first_search(g, vertex(0, g), visitor(OutDegreeVisitor {ns.begin()}));

            // Setup the tree.
            setup_tree(num_nodes, ns.data(), tree_.data());

            auto const qp_in_size = tree_ocp_qp_in_calculate_size(
                num_nodes, nx.data(), nu.data(), nc.data(), tree_.data());

            qp_in_memory_.resize(qp_in_size);
            tree_ocp_qp_in_create(num_nodes, nx.data(), nu.data(), nc.data(), tree_.data(), 
                &qp_in_, qp_in_memory_.data());

            auto const qp_out_size = tree_ocp_qp_out_calculate_size(
                num_nodes, nx.data(), nu.data(), nc.data());

            qp_out_memory_.resize(qp_out_size);
            tree_ocp_qp_out_create(num_nodes, nx.data(), nu.data(), nc.data(),
                &qp_out_, qp_out_memory_.data());

            auto const tdunes_opts_size = treeqp_tdunes_opts_calculate_size(num_nodes);
            tdunes_opts_mem_.resize(tdunes_opts_size);
            treeqp_tdunes_opts_create(num_nodes, &opts_, tdunes_opts_mem_.data());
            treeqp_tdunes_opts_set_default(num_nodes, &opts_);

            for (int ii = 0; ii < num_nodes; ii++)
            {
                opts_.qp_solver[ii] = TREEQP_QPOASES_SOLVER;

                // TODO: in theory, we should set opts->qp_solver[ii] based on the structure of the Hessian
                // matrix for problem ii. But for now we simply use QPOASES because it must always work.
            }

            auto const treeqp_size = treeqp_tdunes_calculate_size(&qp_in_, &opts_);
            qp_solver_memory_.resize(treeqp_size);
            treeqp_tdunes_create(&qp_in_, &opts_, &work_, qp_solver_memory_.data());
        }


        void solve()
        {
            auto const ret = treeqp_tdunes_solve(&qp_in_, &qp_out_, &opts_, &work_);

            if (ret != TREEQP_OPTIMAL_SOLUTION_FOUND)
                throw DualNewtonTreeException(ret);
        }


        void print() const
        {
            tree_ocp_qp_in_print_dims(&qp_in_);
            tree_ocp_qp_in_print(&qp_in_);
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto vertexIndex() const
        {
            return get(boost::vertex_index, graph_);
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
            return iterator_property_map(size_.begin(), vertexIndex());
        }


        auto Q()
        {
            // void tree_ocp_qp_in_set_node_Q_colmajor(const double * const Q, const int lda, tree_ocp_qp_in * const qp_in, const int indx)
            // void tree_ocp_qp_in_get_node_Q_colmajor(double * const Q, const int lda, const tree_ocp_qp_in * const qp_in, const int indx)
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(), size_Q(size()), 
                std::bind(tree_ocp_qp_in_set_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_node_Q_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto R()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(), size_R(size()), 
                std::bind(tree_ocp_qp_in_set_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_node_R_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto S()
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc assumes it to be NX-by-NU.
            // Treating a rowMajor matrix as column-major to treeQP actually transposes it.
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, rowMajor>>(vertexIndex(), size_S(size()), 
                std::bind(tree_ocp_qp_in_set_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_node_S_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto q()
        {
            // void tree_ocp_qp_in_set_node_q(const double * const q, tree_ocp_qp_in * const qp_in, const int indx)
            // void tree_ocp_qp_in_get_node_q(double * const q, const tree_ocp_qp_in * const qp_in, const int indx)
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_x(size()), 
                std::bind(tree_ocp_qp_in_set_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_q, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto r()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_u(size()), 
                std::bind(tree_ocp_qp_in_set_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_r, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto C()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(), size_C(size()), 
                std::bind(tree_ocp_qp_in_set_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_node_C_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto D()
        {
            return detail::makeMatrixPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel, columnMajor>>(vertexIndex(), size_D(size()), 
                std::bind(tree_ocp_qp_in_set_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_node_D_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto ld()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_d(size()), 
                std::bind(tree_ocp_qp_in_set_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_dmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ud()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_d(size()), 
                std::bind(tree_ocp_qp_in_set_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_dmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lx()
        {
            // void tree_ocp_qp_in_set_node_xmin(const double * const xmin, tree_ocp_qp_in * const qp_in, const int indx);
            // void tree_ocp_qp_in_get_node_xmin(double * const xmin, const tree_ocp_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_x(size()), 
                std::bind(tree_ocp_qp_in_set_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_xmin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto ux()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_x(size()), 
                std::bind(tree_ocp_qp_in_set_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_xmax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto lu()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_u(size()), 
                std::bind(tree_ocp_qp_in_set_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_umin, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto uu()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_u(size()), 
                std::bind(tree_ocp_qp_in_set_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_node_umax, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto A()
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_A(size(), graph_), 
                std::bind(tree_ocp_qp_in_set_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_edge_A_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto B()
        {
            return detail::makeMatrixPropertyMap<OcpEdgeDescriptor, DynamicMatrix<Kernel, columnMajor>>(edgeIndex(), size_B(size(), graph_), 
                std::bind(tree_ocp_qp_in_set_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3),
                std::bind(tree_ocp_qp_in_get_edge_B_colmajor, std::placeholders::_1, std::placeholders::_2, &qp_in_, std::placeholders::_3));
        }


        auto b()
        {
            // void tree_ocp_qp_in_set_edge_b(const double * const b, tree_ocp_qp_in * const qp_in, const int indx);
            // void tree_ocp_qp_in_get_edge_b(double * const b, const tree_ocp_qp_in * const qp_in, const int indx);
            return detail::makeVectorPropertyMap<OcpEdgeDescriptor, DynamicVector<Kernel>>(edgeIndex(), size_b(size(), graph_), 
                std::bind(tree_ocp_qp_in_set_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2),
                std::bind(tree_ocp_qp_in_get_edge_b, std::placeholders::_1, &qp_in_, std::placeholders::_2));
        }


        auto x()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_x(size()), 
                std::bind(tree_ocp_qp_out_set_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_ocp_qp_out_get_node_x, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }


        auto u()
        {
            return detail::makeVectorPropertyMap<OcpVertexDescriptor, DynamicVector<Kernel>>(vertexIndex(), size_u(size()), 
                std::bind(tree_ocp_qp_out_set_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2),
                std::bind(tree_ocp_qp_out_get_node_u, std::placeholders::_1, &qp_out_, std::placeholders::_2));
        }

    private:
        OcpGraph graph_;
        std::vector<OcpSize> size_;

        std::vector<node> tree_;
        tree_ocp_qp_in qp_in_;
        std::vector<char> qp_in_memory_;
        tree_ocp_qp_out qp_out_;
        std::vector<char> qp_out_memory_;
        treeqp_tdunes_opts_t opts_;
        std::vector<char> tdunes_opts_mem_;
        treeqp_tdunes_workspace work_;
        std::vector<char> qp_solver_memory_;
    };
}