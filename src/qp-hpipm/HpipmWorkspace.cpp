#include <tmpc/qp/HpipmWorkspace.hpp>

#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <vector>
#include <algorithm>

namespace tmpc
{
	HpipmException::HpipmException(int code)
	:	std::runtime_error("HPIPM return code " + std::to_string(code))
	,	_code(code)
	{
	}
}
