// TODO: get rid of all nasty deprecated stuff

#pragma once

#include "../core/matrix.hpp"

#include <cstddef>
#include <stdexcept>

namespace tmpc
{
	template<typename QP>
	std::size_t constexpr n_xu(QP const& qp) { return qp.nX() + qp.nU(); }

	template<typename QP>
	std::size_t constexpr n_xu() { return QP::nX() + QP::nU(); }

	template<typename QP>
	std::size_t constexpr n_x() { return QP::nX(); }

	template<typename QP>
	std::size_t constexpr n_u() { return QP::nU(); }

	template<typename QP>
	std::size_t constexpr n_d() { return QP::nD(); }

	template<typename QP>
	std::size_t constexpr n_d_end() { return QP::nDT(); }

	template<typename QP>
	typename QP::size_type nIndep(QP const& qp) { return qp.nX() + qp.nU() * qp.nT(); }

	template<typename QP>
	typename QP::size_type nDep(QP const& qp) { return qp.nX() * qp.nT(); }

	template<typename QP>
	typename QP::size_type nVar(QP const& qp) { return qp.nZ() * qp.nT() + qp.nX(); }

	template<typename QP>
	typename QP::size_type nConstr(QP const& qp) { return qp.nD() * qp.nT() + qp.nDT(); }

	template <typename QP>
	Eigen::VectorXd get_xu_min(QP const& qp, std::size_t i)
	{
		auto const x_min = qp.get_x_min(i);
		auto const u_min = qp.get_u_min(i);
		Eigen::VectorXd val(x_min.size() + u_min.size());
		val << x_min, u_min;
		return val;
	}

	template <typename QP, typename Matrix>
	void set_xu_min(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu<QP>(),	"Vector with (NX+NU) rows is expected");
		qp.set_x_min(i, top_rows   <n_x<QP>()>(val));
		qp.set_u_min(i, bottom_rows<n_u<QP>()>(val));
	}

	template<typename QP>
	void set_xu_min(QP& qp, std::size_t i, double val)
	{
		qp.set_x_min(i, QP::StateVector::Constant(val));
		qp.set_u_min(i, QP::InputVector::Constant(val));
	}

	template <typename QP>
	Eigen::VectorXd get_xu_max(QP const& qp, std::size_t i)
	{
		auto const x_max = qp.get_x_max(i);
		auto const u_max = qp.get_u_max(i);
		Eigen::VectorXd val(x_max.size() + u_max.size());
		val << x_max, u_max;
		return val;
	}

	template <typename QP, typename Matrix>
	void set_xu_max(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu<QP>(),	"Vector with (NX+NU) rows is expected");
		qp.set_x_max(i, top_rows   <n_x<QP>()>(val));
		qp.set_u_max(i, bottom_rows<n_u<QP>()>(val));
	}

	template<typename QP>
	void set_xu_max(QP& qp, std::size_t i, double val)
	{
		qp.set_x_max(i, QP::StateVector::Constant(val));
		qp.set_u_max(i, QP::InputVector::Constant(val));
	}

	/**
	 * \brief Set terminal state lower bound
	 */
	template <typename QP, typename Matrix>
	void set_x_end_min(QP& qp, Matrix const& val)
	{
		qp.set_x_min(qp.nT(), val);
	}

	/**
	 * \brief Set terminal state lower bound
	 */
	template <typename QP>
	void set_x_end_min(QP& qp, double val)
	{
		qp.set_x_min(qp.nT(), QP::StateVector::Constant(val));
	}

	/**
	 * \brief Set terminal state upper bound
	 */
	template <typename QP, typename Matrix>
	void set_x_end_max(QP& qp, Matrix const& val)
	{
		qp.set_x_max(qp.nT(), val);
	}

	/**
	 * \brief Set terminal state upper bound
	 */
	template <typename QP>
	void set_x_end_max(QP& qp, double val)
	{
		qp.set_x_max(qp.nT(), QP::StateVector::Constant(val));
	}

	template<typename QP>
	auto get_x_end_min(QP const& qp) -> decltype(qp.get_x_min(qp.nT()))
	{
		return qp.get_x_min(qp.nT());
	}

	template<typename QP>
	auto get_x_end_max(QP const& qp) -> decltype(qp.get_x_max(qp.nT()))
	{
		return qp.get_x_max(qp.nT());
	}

	template <typename QP, typename Matrix>
	void set_H(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu<QP>() && Matrix::ColsAtCompileTime == n_xu<QP>(),
			"Matrix of size (NX+NU)x(NX+NU) is expected");
		qp.set_Q(i, top_left_corner    <n_x<QP>(), n_x<QP>()>(val));
		qp.set_S(i, top_right_corner   <n_x<QP>(), n_u<QP>()>(val));
		qp.set_R(i, bottom_right_corner<n_u<QP>(), n_u<QP>()>(val));
	}

	template <typename QP>
	Eigen::MatrixXd get_H(QP const& qp, std::size_t i)
	{
		auto const Q = qp.get_Q(i);
		auto const R = qp.get_R(i);
		auto const S = qp.get_S(i);
		Eigen::MatrixXd H(Q.rows() + R.rows(), Q.cols() + R.cols());
		H << Q, S, S.transpose(), R;
		return H;
	}

	template <typename QP, typename Matrix>
	void set_g(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix>& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu<QP>(),	"Column vector of size (NX+NU) is expected");
		qp.set_q(i, top_rows   <n_x<QP>()>(val));
		qp.set_r(i, bottom_rows<n_u<QP>()>(val));
	}

	template <typename QP>
	Eigen::VectorXd get_g(QP const& qp, std::size_t i)
	{
		auto const q = qp.get_q(i);
		auto const r = qp.get_r(i);
		Eigen::VectorXd g(q.size() + r.size());
		g << q, r;
		return g;
	}

	template <typename QP, typename Matrix>
	void set_AB(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& AB)
	{
		static_assert(Matrix::ColsAtCompileTime == n_xu<QP>(),	"Matrix with (NX+NU) columns is expected");
		qp.set_A(i, left_cols <n_x<QP>()>(AB));
		qp.set_B(i, right_cols<n_u<QP>()>(AB));
	}

	template <typename QP>
	Eigen::MatrixXd get_AB(QP const& qp, std::size_t i)
	{
		auto const A = qp.get_A(i);
		auto const B = qp.get_B(i);
		Eigen::MatrixXd AB(A.rows(), A.cols() + B.cols());
		AB << A, B;
		return AB;
	}

	template <typename QP, typename Matrix>
	void set_CD(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::ColsAtCompileTime == n_xu<QP>(),	"Matrix with (NX+NU) columns is expected");
		qp.set_C(i, left_cols <n_x<QP>()>(val));
		qp.set_D(i, right_cols<n_u<QP>()>(val));
	}

	template <typename QP>
	Eigen::MatrixXd get_CD(QP const& qp, std::size_t i)
	{
		auto const C = qp.get_C(i);
		auto const D = qp.get_D(i);

		Eigen::MatrixXd CD(C.rows(), C.cols() + D.cols());

		if (CD.rows() > 0)	// workaround for this bug: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1242
			CD << C, D;

		return CD;
	}

	template <typename QP, typename Matrix>
	void set_Q_end(QP& qp, Eigen::MatrixBase<Matrix> const& val)	{
		qp.set_Q(qp.nT(), val);
	}

	template <typename QP>
	auto get_Q_end(QP const& qp) -> decltype(qp.get_Q(qp.nT()))
	{
		return qp.get_Q(qp.nT());
	}

	template <typename QP, typename Matrix>
	void set_q_end(QP& qp, Eigen::MatrixBase<Matrix> const& val)	{
		qp.set_q(qp.nT(), val);
	}

	template <typename QP>
	auto get_q_end(QP const& qp) -> decltype(qp.get_q(qp.nT()))
	{
		return qp.get_q(qp.nT());
	}
}
