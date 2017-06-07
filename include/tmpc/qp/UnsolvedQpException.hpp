/*
 * UnsolvedQpException.hpp
 *
 *  Created on: Dec 5, 2016
 *      Author: mkatliar
 */

#pragma once

#include <stdexcept>

namespace tmpc {

/**
 * \brief Thrown by a QP solver when a QP could not be solved.
 *
 * TODO: Store a QP of generic class instead of a specific type.
 */
class UnsolvedQpException : public std::runtime_error
{
public:
	UnsolvedQpException(std::string const& solver_name) 
	:	std::runtime_error(solver_name + " could not solve a QP.")
	{
	}
};

}	// namespace tmpc
