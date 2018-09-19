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
	}


	DualNewtonTreeException::DualNewtonTreeException(return_t code)
	:	std::runtime_error(dualNewtonTreeExceptionMessage(code))
	,	_code(code)
	{
	}


	DualNewtonTreeOptions::DualNewtonTreeOptions(size_t num_nodes)
	:	mem_(new char[treeqp_tdunes_opts_calculate_size(num_nodes)])
	{
		treeqp_tdunes_opts_create(num_nodes, &opts_, mem_.get());
		treeqp_tdunes_opts_set_default(num_nodes, &opts_);

		for (int ii = 0; ii < num_nodes; ii++)
		{
			opts_.qp_solver[ii] = TREEQP_QPOASES_SOLVER;

			// TODO: in theory, we should set opts->qp_solver[ii] based on the structure of the Hessian
			// matrix for problem ii. But for now we simply use QPOASES because it must always work.
		}
	}
}
