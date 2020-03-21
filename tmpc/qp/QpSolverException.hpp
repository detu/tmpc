/*
 * QpSolverException.hpp
 *
 *  Created on: Dec 5, 2016
 *      Author: mkatliar
 */

#pragma once

#include <tmpc/Exception.hpp>


namespace tmpc 
{
	/**
	 * \brief Thrown by a QP solver when a QP could not be solved.
	 */
	class QpSolverException 
	:	virtual public boost::exception
	,	public std::runtime_error
	{
	public:
		QpSolverException()
		:	std::runtime_error("QP solver could not solve a QP.")
		{
		}
	};
}
