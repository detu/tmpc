#pragma once

#include "QpSize.hpp"

#include <tmpc/Matrix.hpp>

#include <vector>
#include <iostream>

namespace tmpc {

/**
 * \brief Manages input data for qpOASES solver.
 *
 * Implements concept: QuadraticProblem.
 */
class QpOasesProblem
{
	// Matrix storage option -- important!
	// Must be rowMajor, because qpOASES expects input matrices in row-major format.
	static auto constexpr storageOrder = rowMajor;

public:
	typedef unsigned int size_type;
	typedef double Scalar;
	typedef DynamicMatrix<Scalar, storageOrder> Matrix;
	typedef DynamicVector<Scalar, columnVector> Vector;
	typedef Submatrix<Matrix> SubM;
	typedef Subvector<Vector> SubV;

	class Stage
	{
	public:
		Stage(QpSize const& sz, SubM const& Q, SubM const& R, SubM const& S, SubM const& ST, SubV const& q, SubV const& r,
				SubV const& lbx, SubV const& ubx, SubV const& lbu, SubV const& ubu, SubM const& A, SubM const& B, SubV const& lbb, SubV const& ubb,
				SubM const& C, SubM const& D, SubV const& lbd, SubV const& ubd);

		Stage(Stage const&) = delete;
		Stage(Stage &&) = default;

		template <typename Expr>
		Stage& operator=(Expr const& rhs)
		{
			assign(*this, rhs);
			return *this;
		}

		const SubM& get_A() const {
			return A_;
		}

		template <typename T>
		void set_A(const T& a) {
			A_ = a;
		}

		decltype(auto) get_b() const {
			return -lbb_;
		}

		template <typename T>
		void set_b(const T& b) {
			lbb_ = ubb_ = -b;
		}

		const SubM& get_B() const {
			return B_;
		}

		template <typename T>
		void set_B(const T& b) {
			B_ = b;
		}

		const SubM& get_C() const {
			return C_;
		}

		template <typename T>
		void set_C(const T& c) {
			C_ = c;
		}

		const SubM& get_D() const {
			return D_;
		}

		template <typename T>
		void set_D(const T& d) {
			D_ = d;
		}

		const SubV& get_lbd() const {
			return lbd_;
		}

		template <typename T>
		void set_lbd(const T& lbd) {
			lbd_ = lbd;
		}

		const SubV& get_lbu() const {
			return lbu_;
		}

		template <typename T>
		void set_lbu(const T& lbu) {
			lbu_ = lbu;
		}

		const SubV& get_lbx() const {
			return lbx_;
		}

		template <typename T>
		void set_lbx(const T& lbx) {
			lbx_ = lbx;
		}

		const SubV& get_q() const {
			return q_;
		}

		template <typename T>
		void set_q(const T& q) {
			q_ = q;
		}

		const SubM& get_Q() const {
			return Q_;
		}

		template <typename T>
		void set_Q(const T& q) {
			Q_ = q;
		}

		const SubV& get_r() const {
			return r_;
		}

		template <typename T>
		void set_r(const T& r) {
			r_ = r;
		}

		const SubM& get_R() const {
			return R_;
		}

		template <typename T>
		void set_R(const T& r) {
			R_ = r;
		}

		const SubM& get_S() const {
			return S_;
		}

		template <typename T>
		void set_S(const T& s) {
			ST_ = (S_ = s).transpose();
		}

		const SubV& get_ubd() const {
			return ubd_;
		}

		template <typename T>
		void set_ubd(const T& ubd) {
			ubd_ = ubd;
		}

		const SubV& get_ubu() const {
			return ubu_;
		}

		template <typename T>
		void set_ubu(const T& ubu) {
			ubu_ = ubu;
		}

		const SubV& get_ubx() const {
			return ubx_;
		}

		template <typename T>
		void set_ubx(const T& ubx) {
			ubx_ = ubx;
		}

		QpSize const& size() const
		{
			return size_;
		}

	private:
		QpSize size_;
		SubM Q_;
		SubM R_;
		SubM S_;
		SubM ST_;
		SubV q_;
		SubV r_;
		SubV lbx_;
		SubV ubx_;
		SubV lbu_;
		SubV ubu_;
		SubM A_;
		SubM B_;
		SubV lbb_;
		SubV ubb_;
		SubM C_;
		SubM D_;
		SubV lbd_;
		SubV ubd_;
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

	// Is this going to become a "modern" multistage QP interface?
	Stage& operator[](std::size_t i)
	{
		return stage_.at(i);
	}

	Stage const& operator[](std::size_t i) const
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

private:
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

	QpOasesProblem(std::vector<QpSize> const& sz, size_type n_var, size_type n_constr);
	void InitStages();

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
