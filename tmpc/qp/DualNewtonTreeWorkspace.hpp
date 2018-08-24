#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/Matrix.hpp>

#include <treeqp/src/tree_ocp_qp_common.h>
#include <treeqp/src/dual_Newton_tree.h>
#include <treeqp/utils/tree.h>
#include <treeqp/utils/print.h>

#include <boost/graph/breadth_first_search.hpp>

#include <vector>
#include <stdexcept>


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


    template <typename Kernel_>
    class DualNewtonTreeWorkspace
    {
        template <typename Iterator>
        class RecordOcpSizeVisitor
        :   public boost::default_bfs_visitor 
        {
        public:
            RecordOcpSizeVisitor(Iterator nx, Iterator nu, Iterator nc)
            :   nx_{nx}
            ,   nu_{nu}
            ,   nc_{nc}
            {                
            }

        
            template <typename Vertex, typename Graph>
            void discover_vertex(Vertex u, const Graph & g)
            {
                *nx_++ = g[u].size.nx();
                *nu_++ = g[u].size.nu();
                *nc_++ = g[u].size.nc();
            }

        
        private:
            Iterator nx_;
            Iterator nu_;
            Iterator nc_;
        };


        template <typename Iterator>
        static RecordOcpSizeVisitor<Iterator> recordOcpSizeVisitor(Iterator nx, Iterator nu, Iterator nc)
        {
            return RecordOcpSizeVisitor<Iterator>(nx, nu, nc);
        }


    public:
        using Kernel = Kernel_;
        using Real = typename Kernel::Real;
        using vertex_descriptor = boost::graph_traits<OcpSizeGraph>::vertex_descriptor;
        using edge_descriptor = boost::graph_traits<OcpSizeGraph>::edge_descriptor;


        template <
            typename Key,
            std::pair<size_t, size_t> (* SizeFunc)(OcpSizeGraph const&, Key),
            void (* SetterFunc)(const Real * const, const int, tree_ocp_qp_in * const, const int),
            void (* GetterFunc)(Real * const R, const int lda, const tree_ocp_qp_in * const qp_in, const int indx),
            StorageOrder SO = columnMajor
        >
        class MatrixPropertyMap
        {
        public:
            MatrixPropertyMap(DualNewtonTreeWorkspace& ws)
            :   ws_{ws}
            {
            }


            friend void put(MatrixPropertyMap const& pm, Key k, DynamicMatrix<Kernel> const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(MatrixPropertyMap const& pm, Key k)
            {
                return pm.get(k);
            }


        private:
            void put(Key k, DynamicMatrix<Kernel, SO> const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                if (sz.first != rows(val) || sz.second != columns(val))
                    throw std::invalid_argument("Invalid matrix size");

                // TODO: use the proper stride value instead of rows()/columns().
                auto const spacing = SO == columnMajor ? rows(val) : columns(val);
                SetterFunc(val.data(), spacing, &ws_.qp_in_, ws_.getIndex(k));
            }


            DynamicMatrix<Kernel, SO> get(Key k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                
                DynamicMatrix<Kernel, SO> m(sz.first, sz.second);

                // TODO: use the proper stride value instead of rows()/columns().
                auto const spacing = SO == columnMajor ? rows(m) : columns(m);
                GetterFunc(m.data(), spacing, &ws_.qp_in_, ws_.getIndex(k));

                return m;
            }


            DualNewtonTreeWorkspace& ws_;
        };


        template <
            typename Key,
            size_t (* SizeFunc)(OcpSizeGraph const&, Key),
            void (* SetterFunc)(const double * const, tree_ocp_qp_in * const, const int)
        >
        class VectorPropertyMap
        {
        public:
            VectorPropertyMap(DualNewtonTreeWorkspace& ws)
            :   ws_{ws}
            {
            }


            friend void put(VectorPropertyMap const& pm, Key k, DynamicVector<Kernel> const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(VectorPropertyMap const& pm, Key k)
            {
                return pm.get(k);
            }


        private:
            void put(Key k, DynamicVector<Kernel> const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                if (sz != val.size())
                    throw std::invalid_argument("Invalid vector size");

                SetterFunc(val.data(), &ws_.qp_in_, ws_.getIndex(k));
            }


            DynamicVector<Kernel> get(Key k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                // TODO: call getter here
                return DynamicVector<Kernel>(sz);
            }


            DualNewtonTreeWorkspace& ws_;
        };


        template <
            typename Key,
            size_t (* SizeFunc)(OcpSizeGraph const&, Key),
            void (* GetterFunc)(double *, const tree_ocp_qp_out * const, const int)
        >
        class SolutionVectorPropertyMap
        {
        public:
            SolutionVectorPropertyMap(DualNewtonTreeWorkspace& ws)
            :   ws_{ws}
            {
            }


            friend decltype(auto) get(SolutionVectorPropertyMap const& pm, Key k)
            {
                return pm.get(k);
            }


        private:
            DynamicVector<Kernel> get(Key k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                DynamicVector<Kernel> ret(sz);
                GetterFunc(ret.data(), &ws_.qp_out_, ws_.getIndex(k));

                return ret;
            }


            DualNewtonTreeWorkspace& ws_;
        };


        DualNewtonTreeWorkspace(OcpSizeGraph const& g)
        :   graph_{g}
        {
            auto const num_nodes = num_vertices(g);
            tree_.resize(num_nodes);

            // ns: an array of numbers, as if you were traversing the tree breadth-first,
            // and output the number of children for each node.
            // For now we assume linear tree structure.
            std::vector<int> ns(num_nodes - 1, 1);
            ns.push_back(0);
            setup_tree(num_nodes, ns.data(), tree_.data());

            //print_node(tree_.data());

            std::vector<int> nx, nu, nc;
            nx.reserve(num_nodes);
            nu.reserve(num_nodes);
            nc.reserve(num_nodes);
            
            auto vis = recordOcpSizeVisitor(std::back_inserter(nx), std::back_inserter(nu), std::back_inserter(nc));
            breadth_first_search(g, vertex(0, g), visitor(vis));

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


        void setNodeObjective(vertex_descriptor v, 
            DynamicMatrix<Kernel> & Q, DynamicMatrix<Kernel> & R, DynamicMatrix<Kernel> & S,
            DynamicVector<Kernel> & q, DynamicVector<Kernel> & r)
        {
            tree_ocp_qp_in_set_node_objective_colmajor(Q.data(), R.data(), S.data(), q.data(), r.data(), &qp_in_, v);
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
            return get(boost::edge_index, graph_);
        }


        auto size() const
        {
            return get(&OcpVertex::size, std::as_const(graph_));
        }


        auto Q()
        {
            return MatrixPropertyMap<vertex_descriptor, &size_Q, &tree_ocp_qp_in_set_node_Q_colmajor, &tree_ocp_qp_in_get_node_Q_colmajor>(*this);
        }


        auto R()
        {
            return MatrixPropertyMap<vertex_descriptor, &size_R, &tree_ocp_qp_in_set_node_R_colmajor, &tree_ocp_qp_in_get_node_R_colmajor>(*this);
        }


        auto S()
        {
            // NOTE: treeQP assumes the size of S to be NU-by-NX, tmpc assumes it to be NX-by-NU.
            // Passing rowMajor matrix as column-major to treeQP actually transposes it.
            return MatrixPropertyMap<vertex_descriptor, &size_S, &tree_ocp_qp_in_set_node_S_colmajor, &tree_ocp_qp_in_get_node_S_colmajor, rowMajor>(*this);
        }


        auto q()
        {
            return VectorPropertyMap<vertex_descriptor, &size_q, &tree_ocp_qp_in_set_node_q>(*this);
        }


        auto r()
        {
            return VectorPropertyMap<vertex_descriptor, &size_r, &tree_ocp_qp_in_set_node_r>(*this);
        }


        auto C()
        {
            return MatrixPropertyMap<vertex_descriptor, &size_C, &tree_ocp_qp_in_set_node_C_colmajor, &tree_ocp_qp_in_get_node_C_colmajor>(*this);
        }


        auto D()
        {
            return MatrixPropertyMap<vertex_descriptor, &size_D, &tree_ocp_qp_in_set_node_D_colmajor, &tree_ocp_qp_in_get_node_D_colmajor>(*this);
        }


        auto ld()
        {
            return VectorPropertyMap<vertex_descriptor, &size_d, &tree_ocp_qp_in_set_node_dmin>(*this);
        }


        auto ud()
        {
            return VectorPropertyMap<vertex_descriptor, &size_d, &tree_ocp_qp_in_set_node_dmax>(*this);
        }


        auto lx()
        {
            return VectorPropertyMap<vertex_descriptor, &size_x, &tree_ocp_qp_in_set_node_xmin>(*this);
        }


        auto ux()
        {
            return VectorPropertyMap<vertex_descriptor, &size_x, &tree_ocp_qp_in_set_node_xmax>(*this);
        }


        auto lu()
        {
            return VectorPropertyMap<vertex_descriptor, &size_u, &tree_ocp_qp_in_set_node_umin>(*this);
        }


        auto uu()
        {
            return VectorPropertyMap<vertex_descriptor, &size_u, &tree_ocp_qp_in_set_node_umax>(*this);
        }


        auto A()
        {
            return MatrixPropertyMap<edge_descriptor, &size_A, &tree_ocp_qp_in_set_edge_A_colmajor, &tree_ocp_qp_in_get_edge_A_colmajor>(*this);
        }


        auto B()
        {
            return MatrixPropertyMap<edge_descriptor, &size_B, &tree_ocp_qp_in_set_edge_B_colmajor, &tree_ocp_qp_in_get_edge_B_colmajor>(*this);
        }


        auto b()
        {
            return VectorPropertyMap<edge_descriptor, &size_b, &tree_ocp_qp_in_set_edge_b>(*this);
        }


        auto x()
        {
            return SolutionVectorPropertyMap<vertex_descriptor, &size_x, &tree_ocp_qp_out_get_node_x>(*this);
        }


        auto u()
        {
            return SolutionVectorPropertyMap<vertex_descriptor, &size_u, &tree_ocp_qp_out_get_node_u>(*this);
        }


    private:
        //using SizeMap = typename boost::property_map<OcpSizeGraph, OcpSize OcpVertex::*>::const_type;
        int getIndex(vertex_descriptor v) const
        {
            return get(vertexIndex(), v);
        }


        int getIndex(edge_descriptor e) const
        {
            return get(edgeIndex(), e);
        }


        OcpSizeGraph graph_;

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