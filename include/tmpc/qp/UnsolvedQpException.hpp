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

	// Declaring qpData_ as a shared_ptr<> allows copying of UnsolvedQpException.
	// This can be necessary, for example, if we want to start a separate thread to dump a failed QP
	// and at the same time rethrow the exception for further processing.
	// See also https://gitlab.tuebingen.mpg.de/mkatliar/cms/issues/17
	std::shared_ptr<QpData> qpData_;
};

}	// namespace tmpc
