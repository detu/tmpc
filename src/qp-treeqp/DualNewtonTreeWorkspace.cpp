#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>

#include <sstream>
#include <string>
#include <map>


namespace tmpc
{
	namespace
	{
		std::string dualNewtonTreeExceptionMessage(return_t code)
		{
			static const std::map<return_t, char const *> code_to_string =
			{
				{TREEQP_OPTIMAL_SOLUTION_FOUND, "TREEQP_OPTIMAL_SOLUTION_FOUND"},
    			{TREEQP_MAXIMUM_ITERATIONS_REACHED, "TREEQP_MAXIMUM_ITERATIONS_REACHED"},
				{TREEQP_DN_NOT_DESCENT_DIRECTION, "TREEQP_DN_NOT_DESCENT_DIRECTION"},
    			{TREEQP_DN_STAGE_QP_INIT_FAILED, "TREEQP_DN_STAGE_QP_INIT_FAILED"},
				{TREEQP_DN_STAGE_QP_SOLVE_FAILED, "TREEQP_DN_STAGE_QP_SOLVE_FAILED"},
    			{TREEQP_IP_MIN_STEP, "TREEQP_IP_MIN_STEP"},
    			{TREEQP_IP_UNKNOWN_FLAG, "TREEQP_IP_UNKNOWN_FLAG"},
				{TREEQP_OK, "TREEQP_OK"},
				{TREEQP_INVALID_OPTION, "TREEQP_INVALID_OPTION"},
				{TREEQP_ERROR_OPENING_FILE, "TREEQP_ERROR_OPENING_FILE"},
    			{TREEQP_UNKNOWN_ERROR, "TREEQP_UNKNOWN_ERROR"}
			};

			auto const it = code_to_string.find(code);

			std::stringstream msg;
			msg << "treeQP return code " << code << " (" << (it != code_to_string.end() ? it->second : "unrecognized return code") << ")";

			return msg.str();
		}


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
	}


	DualNewtonTreeException::DualNewtonTreeException(return_t code)
	:	std::runtime_error(dualNewtonTreeExceptionMessage(code))
	,	_code(code)
	{
	}


	DualNewtonTreeOptions::DualNewtonTreeOptions(size_t num_nodes)
	#ifdef TMPC_TREEQP_USE_HPMPC
	:	mem_(new char[treeqp_hpmpc_opts_calculate_size(num_nodes)])
	#else
	:	mem_(new char[treeqp_tdunes_opts_calculate_size(num_nodes)])
	#endif
	{
		#ifdef TMPC_TREEQP_USE_HPMPC
		treeqp_hpmpc_opts_create(num_nodes, &opts_, mem_.get());
		treeqp_hpmpc_opts_set_default(num_nodes, &opts_);
		#else
		treeqp_tdunes_opts_create(num_nodes, &opts_, mem_.get());
		treeqp_tdunes_opts_set_default(num_nodes, &opts_);

		for (int ii = 0; ii < num_nodes; ii++)
		{
			opts_.qp_solver[ii] = TREEQP_QPOASES_SOLVER;

			// TODO: in theory, we should set opts->qp_solver[ii] based on the structure of the Hessian
			// matrix for problem ii. But for now we simply use QPOASES because it must always work.
		}
		#endif
	}


	void DualNewtonTreeWorkspace::init()
	{
		auto const num_nodes = num_vertices(graph_);
		auto const vertex_id = vertexIndex(graph_);

		// Fill size arrays.
		std::vector<int> nx(num_nodes), nu(num_nodes), nc(num_nodes), nk(num_nodes);

		for (auto v : vertices(graph_))
		{
			auto const i = vertex_id[v];
			nx[i] = size_[i].nx();
			nu[i] = size_[i].nu();
			nc[i] = size_[i].nc();
		}

		// Fill the number of kids vector with out-degrees of the nodes.
		breadth_first_search(graph_, vertex(0, graph_), visitor(OutDegreeVisitor {nk.begin()}));

		// Setup the tree.
		auto const qp_in_size = tree_qp_in_calculate_size(
			num_nodes, nx.data(), nu.data(), nc.data(), nk.data());

		qp_in_memory_.resize(qp_in_size);
		tree_qp_in_create(num_nodes, nx.data(), nu.data(), nc.data(), nk.data(),
			&qp_in_, qp_in_memory_.data());
	
		// Adding bounds everywhere for HPMPC
		const double inf_val = TREEQP_INF/10.;

		int max_dim = 0;
		for (int ii = 0; ii < num_nodes; ii++)
		{
			if (max_dim < nx[ii]) max_dim = nx[ii];
			if (max_dim < nu[ii]) max_dim = nu[ii];
		}
		double *almost_inf = new double [max_dim];
		double *almost_minus_inf = new double [max_dim];

		for (int ii =0; ii < max_dim; ii++)
		{
			almost_inf[ii] = inf_val;
			almost_minus_inf[ii] = -inf_val;
		}

		for (int ii =0; ii < num_nodes; ii++)
		{
			tree_qp_in_set_node_xmin(almost_minus_inf, &qp_in_, ii);
			tree_qp_in_set_node_xmax(almost_inf, &qp_in_, ii);
			tree_qp_in_set_node_umin(almost_minus_inf, &qp_in_, ii);
			tree_qp_in_set_node_umax(almost_inf, &qp_in_, ii);

		}

		delete[] almost_inf;
		delete[] almost_minus_inf;

		auto const qp_out_size = tree_qp_out_calculate_size(
			num_nodes, nx.data(), nu.data(), nc.data());

		qp_out_memory_.resize(qp_out_size);
		tree_qp_out_create(num_nodes, nx.data(), nu.data(), nc.data(),
			&qp_out_, qp_out_memory_.data());

		#ifdef TMPC_TREEQP_USE_HPMPC
		auto const treeqp_size = treeqp_hpmpc_calculate_size(&qp_in_, &opts_.nativeOptions());
		qp_solver_memory_.resize(treeqp_size);
		treeqp_hpmpc_create(&qp_in_, &opts_.nativeOptions(), &work_, qp_solver_memory_.data());
		#else
		auto const treeqp_size = treeqp_tdunes_calculate_size(&qp_in_, &opts_.nativeOptions());
		qp_solver_memory_.resize(treeqp_size);
		treeqp_tdunes_create(&qp_in_, &opts_.nativeOptions(), &work_, qp_solver_memory_.data());
		#endif
	}


	void DualNewtonTreeWorkspace::print() const
	{
		tree_qp_in_print_dims(&qp_in_);
		tree_qp_in_print(&qp_in_);
	}


	void DualNewtonTreeWorkspace::solve()
	{
		#ifdef TMPC_TREEQP_USE_HPMPC
		auto const ret = treeqp_hpmpc_solve(&qp_in_, &qp_out_, &opts_.nativeOptions(), &work_);
		#else
		// A workaround for this issue: https://gitlab.syscop.de/dimitris.kouzoupis/hangover/issues/6
		if (true)
		{
			size_t sum_nx = 0;
			for (auto v : vertices(graph_))
				sum_nx += get(size(), v).nx();

			std::vector<Real> lambda(sum_nx);
			treeqp_tdunes_set_dual_initialization(lambda.data(), &work_);
		}

		auto const ret = treeqp_tdunes_solve(&qp_in_, &qp_out_, &opts_.nativeOptions(), &work_);
		#endif

		if (ret != TREEQP_OPTIMAL_SOLUTION_FOUND)
			throw DualNewtonTreeException(ret);
	}
}
