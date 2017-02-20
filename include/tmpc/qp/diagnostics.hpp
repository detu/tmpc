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

	for (unsigned i = 0; i < qp.nT(); ++i)
	{
		CheckFinite("Q", qp.get_Q(i), i, msg_out);
		CheckFinite("R", qp.get_R(i), i, msg_out);
		CheckFinite("S", qp.get_S(i), i, msg_out);
		CheckFinite("q", qp.get_q(i), i, msg_out);
		CheckFinite("r", qp.get_r(i), i, msg_out);
		CheckFinite("A", qp.get_A(i), i, msg_out);
		CheckFinite("B", qp.get_B(i), i, msg_out);
		CheckFinite("C", qp.get_C(i), i, msg_out);
		CheckFinite("D", qp.get_D(i), i, msg_out);
		CheckNotNaN("x_min", qp.get_x_min(i), i, msg_out);
		CheckNotNaN("x_max", qp.get_x_max(i), i, msg_out);
		CheckNotNaN("u_min", qp.get_u_min(i), i, msg_out);
		CheckNotNaN("u_max", qp.get_u_max(i), i, msg_out);
		CheckNotNaN("d_min", qp.get_d_min(i), i, msg_out);
		CheckNotNaN("d_max", qp.get_d_max(i), i, msg_out);
	}

	CheckFinite("Q_end", get_Q_end(qp), qp.nT(), msg_out);
	CheckFinite("q_end", get_q_end(qp), qp.nT(), msg_out);
	CheckFinite("C_end", qp.get_C_end(), qp.nT(), msg_out);
	CheckNotNaN("x_end_min", get_x_end_min(qp), qp.nT(), msg_out);
	CheckNotNaN("x_end_max", get_x_end_max(qp), qp.nT(), msg_out);
	CheckNotNaN("d_end_min", qp.get_d_end_min(), qp.nT(), msg_out);
	CheckNotNaN("d_end_max", qp.get_d_end_max(), qp.nT(), msg_out);

	return msg_out;
}

}	// namespace qp
}	// namespace diagnostics
