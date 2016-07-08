#pragma once

#include <Eigen/Dense>

#include <cstddef>

namespace tmpc
{
	class MultiStageQPSolutionTag {};
	class MultiStageQPTag {};

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
	Eigen::Matrix<double, n_xu<QP>(), 1> get_xu_min(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_xu(qp), 1> val;
		val << qp.get_x_min(i), qp.get_u_min(i);
		return val;
	}

	template <typename QP, typename Matrix>
	void set_xu_min(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu(qp),	"Vector with (NX+NU) rows is expected");
		qp.set_x_min(i, val.template topRows   <qp.nX()>());
		qp.set_u_min(i, val.template bottomRows<qp.nU()>());
	}

	template<typename QP>
	void set_xu_min(QP& qp, std::size_t i, double val)
	{
		qp.set_x_min(i, QP::StateVector::Constant(val));
		qp.set_u_min(i, QP::InputVector::Constant(val));
	}

	template <typename QP>
	Eigen::Matrix<double, n_xu<QP>(), 1> get_xu_max(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_xu(qp), 1> val;
		val << qp.get_x_max(i), qp.get_u_max(i);
		return val;
	}

	template <typename QP, typename Matrix>
	void set_xu_max(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu(qp),	"Vector with (NX+NU) rows is expected");
		qp.set_x_max(i, val.template topRows   <qp.nX()>());
		qp.set_u_max(i, val.template bottomRows<qp.nU()>());
	}

	template<typename QP>
	void set_xu_max(QP& qp, std::size_t i, double val)
	{
		qp.set_x_max(i, QP::StateVector::Constant(val));
		qp.set_u_max(i, QP::InputVector::Constant(val));
	}

	template<typename QP>
	void set_x_end_min(QP& qp, double val)
	{
		qp.set_x_min(qp.nT(), QP::StateVector::Constant(val));
	}

	template<typename QP>
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
		static_assert(Matrix::RowsAtCompileTime == n_xu(qp) && Matrix::ColsAtCompileTime == n_xu(qp),
			"Matrix of size (NX+NU)x(NX+NU) is expected");
		qp.set_Q(i, val.template topLeftCorner    <qp.nX(), qp.nX()>());
		qp.set_S(i, val.template topRightCorner   <qp.nX(), qp.nU()>());
		qp.set_R(i, val.template bottomRightCorner<qp.nU(), qp.nU()>());
	}

	template <typename QP>
	Eigen::Matrix<double, n_xu<QP>(), n_xu<QP>()> get_H(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_xu(qp), n_xu(qp)> H;
		H << qp.get_Q(i), qp.get_S(i), qp.get_S(i).transpose(), qp.get_R(i);
		return H;
	}

	template <typename QP, typename Matrix>
	void set_g(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix>& val)
	{
		static_assert(Matrix::RowsAtCompileTime == n_xu(qp),	"Column vector of size (NX+NU) is expected");
		qp.set_q(i, val.template topRows   <qp.nX()>());
		qp.set_r(i, val.template bottomRows<qp.nU()>());
	}

	template <typename QP>
	Eigen::Matrix<double, n_xu<QP>(), 1> get_g(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_xu(qp), 1> g;
		g << qp.get_q(i), qp.get_r(i);
		return g;
	}

	template <typename QP, typename Matrix>
	void set_AB(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& AB)
	{
		static_assert(Matrix::ColsAtCompileTime == n_xu(qp),	"Matrix with (NX+NU) columns is expected");
		qp.set_A(i, AB.template leftCols <qp.nX()>());
		qp.set_B(i, AB.template rightCols<qp.nU()>());
	}

	template <typename QP>
	Eigen::Matrix<double, n_x<QP>(), n_xu<QP>()> get_AB(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_x<QP>(), n_xu<QP>()> AB;
		AB << qp.get_A(i), qp.get_B(i);
		return AB;
	}

	template <typename QP, typename Matrix>
	void set_CD(QP& qp, std::size_t i, Eigen::MatrixBase<Matrix> const& val)
	{
		static_assert(Matrix::ColsAtCompileTime == n_xu(qp),	"Matrix with (NX+NU) columns is expected");
		qp.set_C(i, val.template leftCols <qp.nX()>());
		qp.set_D(i, val.template rightCols<qp.nU()>());
	}

	template <typename QP>
	Eigen::Matrix<double, n_d<QP>(), n_xu<QP>()> get_CD(QP const& qp, std::size_t i)
	{
		Eigen::Matrix<double, n_d<QP>(), n_xu<QP>()> CD;

		if (qp.nD() > 0)	// workaround for this bug: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1242
			CD << qp.get_C(i), qp.get_D(i);

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
