#include <qp/CondensingSolver.hpp>

#include <sstream>
#include <stdexcept>

namespace tmpc {
namespace detail {

qpOASES::Options qpOASES_DefaultOptions()
{
	qpOASES::Options options;
	options.setToMPC();
	options.printLevel = qpOASES::PL_LOW;

	assert(options.ensureConsistency() == qpOASES::SUCCESSFUL_RETURN);

	return options;
}

}	// namespace detail
}	// namespace tmpc
