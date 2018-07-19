#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/Matrix.hpp>

#include <treeqp/src/tree_ocp_qp_common.h>
#include <treeqp/src/dual_Newton_tree.h>
#include <treeqp/utils/tree.h>

#include <boost/graph/breadth_first_search.hpp>

#include <vector>


namespace tmpc
{
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
        using vertex_descriptor = boost::graph_traits<OcpSizeGraph>::vertex_descriptor;


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
            treeqp_tdunes_solve(&qp_in_, &qp_out_, &opts_, &work_);
        }


        void print() const
        {
            //tree_ocp_qp_in_print(const_cast<tree_ocp_qp_in *>(&qp_in_));
        }


    private:
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