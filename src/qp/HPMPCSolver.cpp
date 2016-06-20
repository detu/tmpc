/*
 * HPMPCSolver.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: kotlyar
 */

#include <stdexcept>
#include <sstream>

namespace tmpc
{
	void throw_hpmpc_error(int err_code)
	{
		std::ostringstream msg;
		msg << "HPMPC error: return code = " << err_code;
		throw std::runtime_error(msg.str());
	}
}
