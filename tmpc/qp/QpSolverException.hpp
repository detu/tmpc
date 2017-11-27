/*
 * QpSolverException.hpp
 *
 *  Created on: Dec 5, 2016
 *      Author: mkatliar
 */

#pragma once

#include <stdexcept>


namespace tmpc 
{
	/**
	 * \brief Thrown by a QP solver when a QP could not be solved.
	 */
	class QpSolverException : public std::runtime_error
	{
	public:
		QpSolverException(std::string const& solver_name) 
		:	std::runtime_error(solver_name + " could not solve a QP.")
		,	solverName_(solver_name)
		{
		}

		
		QpSolverException(std::string const& solver_name, std::string const& message) 
		:	std::runtime_error(message)
		,	solverName_(solver_name)
		{
		}


		std::string const& solverName() const
		{
			return solverName_;
		}

	private:
		std::string const solverName_;
	};
}
