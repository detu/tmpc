#pragma once

#include <tmpc/Matrix.hpp>

#include <sstream>
#include <cmath>

namespace tmpc {
namespace qp {

namespace detail {

template <typename MT, bool SO>
bool isfinite(Matrix<MT, SO> const& m)
{
	for (size_t i = 0; i < SO == rowMajor ? rows(m) : columns(m); ++i)
	{
		if (std::find_if(m.begin(i), m.end(i),
			[] (typename MT::ElementType const& x) { return std::isfinite(x); }) != m.end(i))
			return false;
	}

	return true;
}

template <typename Matrix, typename OutputIterator>
void CheckFinite(std::string const& name, Matrix const& m, unsigned stage_index, OutputIterator& msg_out)
{
	if (!isfinite(m))
	{
		std::stringstream msg;
		msg << "On stage " << stage_index << ", " << name << " is not all finite:" << std::endl << m;
		*msg_out++ = msg.str();
	}
}

template <typename Matrix, typename OutputIterator>
void CheckNotNaN(std::string const& name, Matrix const& m, unsigned stage_index, OutputIterator& msg_out)
{
	if (isnan(m))
	{
		std::stringstream msg;
		msg << "On stage " << stage_index << ", " << name << " has NaN:" << std::endl << m;
		*msg_out++ = msg.str();
	}
}

}	// namespace detail

/**
 * \brief Checks if a QP is well-formed.
 *
 * \tparam QP class implementing a MultistageQP concept
 * \tparam MessageOutputIterator output iterator for string messages
 *
 * \param message_iterator an iterator to output messages. *message_iterator++ = std::string("string"); should be a valid expression.
 *
 * There should be no NaNs in the data. Infs are allowed only in the bounds.
 *
 * TODO: add lb <= ub checking.
 * TODO: check positive-definiteness of the matrices.\
 * TODO: check feasibility of the constraints?
 */
template <typename QP, typename OutputIterator>
OutputIterator diagnose(QP const& qp, OutputIterator msg_out)
{
	using detail::CheckFinite;
	using detail::CheckNotNaN;

	size_t i = 0;
	for (auto const& st : qp)
	{
		CheckFinite("Q", st.Q(), i, msg_out);
		CheckFinite("R", st.R(), i, msg_out);
		CheckFinite("S", st.S(), i, msg_out);
		CheckFinite("q", st.q(), i, msg_out);
		CheckFinite("r", st.r(), i, msg_out);
		CheckFinite("A", st.A(), i, msg_out);
		CheckFinite("B", st.B(), i, msg_out);
		CheckFinite("C", st.C(), i, msg_out);
		CheckFinite("D", st.D(), i, msg_out);
		CheckNotNaN("x_min", st.lbx(), i, msg_out);
		CheckNotNaN("x_max", st.ubx(), i, msg_out);
		CheckNotNaN("u_min", st.lbu(), i, msg_out);
		CheckNotNaN("u_max", st.ubu(), i, msg_out);
		CheckNotNaN("d_min", st.lbd(), i, msg_out);
		CheckNotNaN("d_max", st.ubd(), i, msg_out);

		++i;
	}

	return msg_out;
}

}	// namespace qp
}	// namespace diagnostics
