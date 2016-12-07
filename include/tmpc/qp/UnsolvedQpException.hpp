/*
 * UnsolvedQpException.hpp
 *
 *  Created on: Dec 5, 2016
 *      Author: mkatliar
 */

#pragma once

#include "Printing.hpp"

#include <memory>

namespace tmpc {

/**
 * \brief Thrown by a QP solver when a QP could not be solved.
 *
 * TODO: Store a QP of generic class instead of a specific type.
 */
class UnsolvedQpException : public std::runtime_error
{
public:
	template <typename QP>
	UnsolvedQpException(std::string const& solver_name, QP const& qp) :
		std::runtime_error(solver_name + " could not solve a QP."),
		qpData_(new QpDataImpl<QP>(qp))
	{
	}

	void PrintQpAsMatlab(std::ostream& os) const
	{
		qpData_->PrintQpAsMatlab(os);
	}

private:
	class QpData
	{
	public:
		virtual ~QpData() {}
		virtual void PrintQpAsMatlab(std::ostream& os) const = 0;
	};

	template <typename QP>
	class QpDataImpl : public QpData
	{
	public:
		QpDataImpl(QP const& qp) : qp_(qp) {}

		void PrintQpAsMatlab(std::ostream& os) const override
		{
			PrintMultistageQpMatlab(os, qp_, "qp");
		}

	private:
		/// \brief Original multistage QP that failed.
		QP const qp_;
	};

	std::unique_ptr<QpData> qpData_;
};

}	// namespace tmpc
