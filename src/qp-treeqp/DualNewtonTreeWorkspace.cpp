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
}
