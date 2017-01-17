#pragma once

#include "QpSize.hpp"

#include <Eigen/Dense>

#include <initializer_list>
#include <vector>
#include <iostream>

namespace tmpc {

/**
 * \brief Manages input data for qpOASES solver.
 *
 * Implements concept: QuadraticProgram.
 */
class QpOasesProblem
{
	// Matrix storage option for Eigen -- important!
	// Must be RowMajor, because qpOASES expects input matrices in column-major format.
	static const int Options = Eigen::RowMajor;

public:
	typedef unsigned int size_type;
	typedef double Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
	typedef Eigen::Map<Matrix, Eigen::Unaligned, Eigen::OuterStride<>> MatrixMap;
	typedef Eigen::Map<Vector> VectorMap;

	class Stage
	{
	public:
		Stage(QpSize const& sz, std::size_t nx_next, std::size_t stride,
				Scalar * H, Scalar * g,	Scalar * lb, Scalar * ub, Scalar * A, Scalar * lbA, Scalar * ubA);

		Stage(Stage const&) = delete;
		Stage(Stage &&) = default;

		const MatrixMap& get_A() const {
			return A_;
		}

		template <typename T>
		void set_A(const Eigen::MatrixBase<T>& a) {
			A_ = a;
		}

		decltype(auto) get_b() const {
			return -lbb_;
		}

		template <typename T>
		void set_b(const Eigen::MatrixBase<T>& b) {
			lbb_ = ubb_ = -b;
		}

		const MatrixMap& get_B() const {
			return B_;
		}

		template <typename T>
		void set_B(const Eigen::MatrixBase<T>& b) {
			B_ = b;
		}

		const MatrixMap& get_C() const {
			return C_;
		}

		template <typename T>
		void set_C(const Eigen::MatrixBase<T>& c) {
			C_ = c;
		}

		const MatrixMap& get_D() const {
			return D_;
		}

		template <typename T>
		void set_D(const Eigen::MatrixBase<T>& d) {
			D_ = d;
		}

		const VectorMap& get_lbd() const {
			return lbd_;
		}

		template <typename T>
		void set_lbd(const Eigen::MatrixBase<T>& lbd) {
			lbd_ = lbd;
		}

		const VectorMap& get_lbu() const {
			return lbu_;
		}

		template <typename T>
		void set_lbu(const Eigen::MatrixBase<T>& lbu) {
			lbu_ = lbu;
		}

		const VectorMap& get_lbx() const {
			return lbx_;
		}

		template <typename T>
		void set_lbx(const Eigen::MatrixBase<T>& lbx) {
			lbx_ = lbx;
		}

		const VectorMap& get_q() const {
			return q_;
		}

		template <typename T>
		void set_q(const Eigen::MatrixBase<T>& q) {
			q_ = q;
		}

		const MatrixMap& get_Q() const {
			return Q_;
		}

		template <typename T>
		void set_Q(const Eigen::MatrixBase<T>& q) {
			Q_ = q;
		}

		const VectorMap& get_r() const {
			return r_;
		}

		template <typename T>
		void set_r(const Eigen::MatrixBase<T>& r) {
			r_ = r;
		}

		const MatrixMap& get_R() const {
			return R_;
		}

		template <typename T>
		void set_R(const Eigen::MatrixBase<T>& r) {
			R_ = r;
		}

		const MatrixMap& get_S() const {
			return S_;
		}

		template <typename T>
		void set_S(const Eigen::MatrixBase<T>& s) {
			ST_ = (S_ = s).transpose();
		}

		const VectorMap& get_ubd() const {
			return ubd_;
		}

		template <typename T>
		void set_ubd(const Eigen::MatrixBase<T>& ubd) {
			ubd_ = ubd;
		}

		const VectorMap& get_ubu() const {
			return ubu_;
		}

		template <typename T>
		void set_ubu(const Eigen::MatrixBase<T>& ubu) {
			ubu_ = ubu;
		}

		const VectorMap& get_ubx() const {
			return ubx_;
		}

		template <typename T>
		void set_ubx(const Eigen::MatrixBase<T>& ubx) {
			ubx_ = ubx;
		}

		QpSize const& size() const
		{
			return size_;
		}

	private:
		QpSize size_;
		MatrixMap Q_;
		MatrixMap R_;
		MatrixMap S_;
		MatrixMap ST_;
		VectorMap q_;
		VectorMap r_;
		VectorMap lbx_;
		VectorMap ubx_;
		VectorMap lbu_;
		VectorMap ubu_;
		MatrixMap A_;
		MatrixMap B_;
		VectorMap lbb_;
		VectorMap ubb_;
		MatrixMap C_;
		MatrixMap D_;
		VectorMap lbd_;
		VectorMap ubd_;
	};

	QpOasesProblem(size_type nx, size_type nc);
	/**
	 * \brief Initialize from QpSize initializer list.
	 */
	QpOasesProblem(std::initializer_list<QpSize> sz);

	/**
	 * \brief Initialize from a vector of QpSize.
	 */
	QpOasesProblem(std::vector<QpSize> const& sz);

	/**
	 * \brief Initialize from QpSize list defined by iterator range.
	 */
	template <typename InputIterator>
	QpOasesProblem(InputIterator size_begin, InputIterator size_end)
	:	QpOasesProblem(std::vector<QpSize>(size_begin, size_end))
	{
	}

	QpOasesProblem(QpOasesProblem const& rhs);
	QpOasesProblem(QpOasesProblem &&) = default;

	size_type nx() const { return static_cast<size_type>(_H.rows()); }
	size_type nc() const { return static_cast<size_type>(_A.rows()); }
	size_type nT() const { return nT_; }

	// Is this going to become a "modern" multistage QP interface?
	Stage& operator[](std::size_t i)
	{
		return stage_.at(i);
	}

	Stage const& operator[](std::size_t i) const
	{
		return stage_.at(i);
	}

	Stage& stage(std::size_t i)
	{
		return stage_.at(i);
	}

	Stage const& stage(std::size_t i) const
	{
		return stage_.at(i);
	}

	std::size_t size() const
	{
		return stage_.size();
	}

	typedef std::vector<Stage>::iterator iterator;
	typedef std::vector<Stage>::const_iterator const_iterator;
	typedef std::vector<Stage>::reference reference;
	typedef std::vector<Stage>::const_reference const_reference;

	iterator begin()
	{
		return stage_.begin();
	}

	iterator end()
	{
		return stage_.end();
	}

	const_iterator begin() const
	{
		return stage_.begin();
	}

	const_iterator end() const
	{
		return stage_.end();
	}

	reference front()
	{
		return stage_.front();
	}

	reference back()
	{
		return stage_.back();
	}

	const_reference front() const
	{
		return stage_.front();
	}

	const_reference back() const
	{
		return stage_.back();
	}

	// TODO: declare deprecated?
	std::vector<QpSize> const& stageSize() const
	{
		return size_;
	}

	// ******************************************************
	//
	// Multistage QP interface
	//
	// ******************************************************
	MatrixMap const& get_Q(size_type i) const
	{
		return stage(i).get_Q();
	}

	template <typename Matrix>
	void set_Q(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_Q(val);
	}

	/**
	 * \brief Get R matrix of stage k
	 */
	MatrixMap const& get_R(size_type k) const
	{
		return stage(k).get_R();
	}

	/**
	 * \brief Set R matrix of a given stage
	 */
	template <typename Matrix>
	void set_R(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).set_R(val);
	}

	/**
	 * \brief Get S matrix of stage k
	 */
	MatrixMap const& get_S(size_type k) const
	{
		return stage(k).get_S();
	}

	/**
	 * \brief Set S matrix of stage k
	 */
	template <typename Matrix>
	void set_S(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).set_S(val);
	}

	/**
	 * \brief Get q vector of stage i
	 */
	VectorMap const& get_q(size_type i) const
	{
		return stage(i).get_q();
	}

	/**
	 * \brief Set q vector of stage i
	 */
	template <typename Matrix>
	void set_q(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_q(val);
	}

	/**
	 * \brief Get r vector of stage k
	 */
	VectorMap const& get_r(size_type k) const
	{
		return stage(k).get_r();
	}

	/**
	 * \brief Set r vector of stage k
	 */
	template <typename Matrix>
	void set_r(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).set_r(val);
	}

	/**
	 * \brief Get A matrix of stage i
	 */
	MatrixMap const& get_A(size_type i) const
	{
		return stage(i).get_A();
	}

	/**
	 * \brief Set A matrix of stage i
	 */
	template <typename Matrix>
	void set_A(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_A(val);
	}

	/**
	 * \brief Get B matrix of stage k
	 */
	MatrixMap const& get_B(size_type k) const
	{
		return stage(k).get_B();
	}

	/**
	 * \brief Set B matrix of stage k
	 */
	template <typename Matrix>
	void set_B(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).set_B(val);
	}

	decltype(auto) get_b(size_type i) const
	{
		return stage(i).get_b();
	}

	template <typename Matrix>
	void set_b(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_b(val);
	}

	MatrixMap const& get_C(size_type i) const
	{
		return stage(i).get_C();
	}

	template <typename Matrix>
	void set_C(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_C(val);
	}

	MatrixMap const& get_D(size_type i) const
	{
		return stage(i).get_D();
	}

	template <typename Matrix>
	void set_D(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_D(val);
	}

	MatrixMap const& get_C_end() const
	{
		return stage_.back().get_C();
	}

	template <typename Matrix>
	void set_C_end(Eigen::MatrixBase<Matrix> const& val)
	{
		stage_.back().set_C(val);
	}

	VectorMap const& get_d_min(size_type i) const
	{
		return stage(i).get_lbd();
	}

	template <typename Matrix>
	void set_d_min(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_lbd(val);
	}

	VectorMap const& get_d_end_min() const
	{
		return stage_.back().get_lbd();
	}

	template <typename Matrix>
	void set_d_end_min(Eigen::MatrixBase<Matrix> const& val)
	{
		stage_.back().set_lbd(val);
	}

	VectorMap const& get_d_max(size_type i) const
	{
		return stage(i).get_ubd();
	}

	template <typename Matrix>
	void set_d_max(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_ubd(val);
	}

	VectorMap const& get_d_end_max() const
	{
		return stage_.back().get_ubd();
	}

	template <typename Matrix>
	void set_d_end_max(Eigen::MatrixBase<Matrix> const& val)
	{
		stage_.back().set_ubd(val);
	}

	VectorMap const& get_x_min(size_type i) const
	{
		return stage(i).get_lbx();
	}

	template <typename Matrix>
	void set_x_min(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_lbx(val);
	}

	VectorMap const& get_x_max(size_type i) const
	{
		return stage(i).get_ubx();
	}

	template <typename Matrix> void set_x_max(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_ubx(val);
	}

	VectorMap const& get_u_min(size_type i) const
	{
		return stage(i).get_lbu();
	}

	template <typename Matrix>
	void set_u_min(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_lbu(val);
	}

	VectorMap const& get_u_max(size_type i) const
	{
		return stage(i).get_ubu();
	}

	template <typename Matrix>
	void set_u_max(size_type i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(i).set_ubu(val);
	}

	//
	// Full matrix and vector access functions.
	//
	Matrix& H() { return _H; }
	const Matrix& H() const { return _H; }

	Vector& g() { return _g; }
	const Vector& g() const { return _g; }

	Matrix& A() { return _A; }
	const Matrix& A() const { return _A; }

	Vector& lbA() { return _lbA; }
	const Vector& lbA() const { return _lbA; }

	Vector& ubA() { return _ubA; }
	const Vector& ubA() const { return _ubA; }

	Vector& lb() { return _lb; }
	const Vector& lb() const { return _lb; }

	Vector& ub() { return _ub; }
	const Vector& ub() const { return _ub; }

	//
	// Raw data access functions for qpOASES.
	//
	double const * H_data() const noexcept { return _H.data(); }
	double const * g_data() const noexcept { return _g.data(); }
	double const * A_data() const noexcept { return _A.data(); }
	double const * lbA_data() const noexcept { return _lbA.data(); }
	double const * ubA_data() const noexcept { return _ubA.data(); }
	double const * lb_data() const noexcept { return _lb.data(); }
	double const * ub_data() const noexcept { return _ub.data(); }

private:
	QpOasesProblem(std::vector<QpSize> const& sz, size_type n_var, size_type n_constr);
	void InitStages();

	std::vector<QpSize> const size_;
	size_type nT_;
	std::vector<Stage> stage_;

	Matrix _H;
	Vector _g;

	// The layout of _lb is [lbx, lbu, ...]
	Vector _lb;

	// The layout of _ub is [ubx, ubu, ...]
	Vector _ub;

	Matrix _A;

	// The layout of _lbA is [lbb, lbd, ...]
	Vector _lbA;

	// The layout of _ubA is [ubb, ubd, ...]
	Vector _ubA;
};

void Print_MATLAB(std::ostream& log_stream, QpOasesProblem const& qp, std::string const& var_name);

}	// namespace tmpc
