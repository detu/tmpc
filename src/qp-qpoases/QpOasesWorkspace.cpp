#include <tmpc/qp/QpOasesWorkspace.hpp>

#include <sstream>
#include <stdexcept>
#include <limits>

namespace tmpc 
{
	static qpOASES::Options qpOASES_DefaultOptions()
	{
		qpOASES::Options options;
		options.setToMPC();
		options.printLevel = qpOASES::PL_LOW;

		assert(options.ensureConsistency() == qpOASES::SUCCESSFUL_RETURN);

		return options;
	}

	namespace detail 
	{
		qpOASES::Options qpOASES_DefaultOptions()
		{
			return ::tmpc::qpOASES_DefaultOptions();
		}
	}

	QpOasesException::QpOasesException(qpOASES::returnValue code)
	:	std::runtime_error("qpOASES return code " + std::to_string(code))
	,	_code(code)
	{
	}
}
