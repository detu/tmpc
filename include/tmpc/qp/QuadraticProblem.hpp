#pragma once

#include <tmpc/Matrix.hpp>

#include "QpSize.hpp"
#include "../core/RealtimeIteration.hpp"	// for RtiQpSize only

#include <vector>
#include <initializer_list>

namespace tmpc {

/** \brief Stores data for a multistage QP problem.
 *
 * Storage format is not explicitly defined and no access to raw data is provided..
 *
 *	TODO: Update problem statement.
 *
 *  The problem is stated as following:
*
*	min  sum_{ k = 0..nI } z_k'*H_k*z_k + g_k'*z_k
*	s.t. x_{ k + 1 } = C_k * z_k + c_k				for k = 0..nI - 1
*            dLow_k <= D_k * z_k <= dUpp_k			for k = 0..nI
*            zMin_k <= z_k <= zMax_k                for k = 0..nI
*
*	where x_k is implicitly defined by z_k = [x_k  u_k] as the first nX variables of z_k
*
*	It holds
*	z_k  \in R^nZ  for k = 0..nI - 1
*   z_nI \in R*nX
*
*	nX < nZ
*	nU = nZ - nX
*/
template <typename Scalar_>
class QuadraticProblem
{
public:
	typedef std::size_t size_type;
	typedef Scalar_ Scalar;
	typedef DynamicMatrix<Scalar> Matrix;
	typedef DynamicVector<Scalar> Vector;

	class Stage
	{
	public:
		Stage(QpSize const& sz, size_type nx_next)
		:	size_(sz)
		,	Q_(sz.nx(), sz.nx())
		,	R_(sz.nu(), sz.nu())
		,	S_(sz.nx(), sz.nu())
		,	q_(sz.nx())
		,	r_(sz.nu())
		,	A_(nx_next, sz.nx())
		,	B_(nx_next, sz.nu())
		,	b_(nx_next)
		,	C_(sz.nc(), sz.nx())
		,	D_(sz.nc(), sz.nu())
		,	lbx_(sz.nx())
		,	ubx_(sz.nx())
		,	lbu_(sz.nu())
		,	ubu_(sz.nu())
		,	lbd_(sz.nc())
		,	ubd_(sz.nc())
		{
		}

		Stage(Stage const&) = delete;
		Stage(Stage &&) = default;

		template <typename Expr>
		Stage& operator=(Expr const& rhs)
		{
			assign(*this, rhs);
			return *this;
		}

		const Matrix& get_A() const {
			return A_;
		}

		template <typename T>
		void set_A(const T& a) {
			A_ = a;
		}

		const Matrix& get_B() const {
			return B_;
		}

		template <typename T>
		void set_B(const T& b) {
			B_ = b;
		}

		Vector const& get_b() const {
			return b_;
		}

		template <typename T>
		void set_b(const T& b) {
			b_ = b;
		}

		const Matrix& get_C() const {
			return C_;
		}

		template <typename T>
		void set_C(const T& c) {
			C_ = c;
		}

		const Matrix& get_D() const {
			return D_;
		}

		template <typename T>
		void set_D(const T& d) {
			D_ = d;
		}

		const Vector& get_lbd() const {
			return lbd_;
		}

		template <typename T>
		void set_lbd(const T& lbd) {
			lbd_ = lbd;
		}

		const Vector& get_lbu() const {
			return lbu_;
		}

		template <typename T>
		void set_lbu(const T& lbu) {
			lbu_ = lbu;
		}

		const Vector& get_lbx() const {
			return lbx_;
		}

		template <typename T>
		void set_lbx(const T& lbx) {
			lbx_ = lbx;
		}

		const Matrix& get_Q() const {
			return Q_;
		}

		template <typename T>
		void set_Q(const T& q) {
			Q_ = q;
		}

		const Matrix& get_R() const {
			return R_;
		}

		template <typename T>
		void set_R(const T& r) {
			R_ = r;
		}

		const Matrix& get_S() const {
			return S_;
		}

		template <typename T>
		void set_S(const T& s) {
			S_ = s;
		}

		const Vector& get_q() const {
			return q_;
		}

		template <typename T>
		void set_q(const T& q) {
			q_ = q;
		}

		const Vector& get_r() const {
			return r_;
		}

		template <typename T>
		void set_r(const T& r) {
			r_ = r;
		}

		const Vector& get_ubd() const {
			return ubd_;
		}

		template <typename T>
		void set_ubd(const T& ubd) {
			ubd_ = ubd;
		}

		const Vector& get_ubu() const {
			return ubu_;
		}

		template <typename T>
		void set_ubu(const T& ubu) {
			ubu_ = ubu;
		}

		const Vector& get_ubx() const {
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
		Matrix Q_;
		Matrix R_;
		Matrix S_;
		Vector q_;
		Vector r_;
		Matrix A_;
		Matrix B_;
		Vector b_;
		Matrix C_;
		Matrix D_;
		Vector lbx_;
		Vector ubx_;
		Vector lbu_;
		Vector ubu_;
		Vector lbd_;
		Vector ubd_;
	};

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

	typedef typename std::vector<Stage>::iterator iterator;
	typedef typename std::vector<Stage>::const_iterator const_iterator;
	typedef typename std::vector<Stage>::reference reference;
	typedef typename std::vector<Stage>::const_reference const_reference;

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

	QuadraticProblem(std::initializer_list<QpSize> sz)
	:	QuadraticProblem(sz.begin(), sz.end())
	{
	}

	template <typename InIter>
	QuadraticProblem(InIter sz_begin, InIter sz_end)
	{
		stage_.reserve(std::distance(sz_begin, sz_end));

		for (auto sz = sz_begin; sz != sz_end; ++sz)
		{
			auto const nx_next = sz + 1 != sz_end ? (sz + 1)->nx() : 0;
			stage_.emplace_back(*sz, nx_next);
		}
	}

	// Default copy constructor is ok.
	QuadraticProblem(QuadraticProblem const&) = default;

private:
	// Private data members.
	//

	// Stores stage data
	std::vector<Stage> stage_;
};

}	// namespace tmpc
