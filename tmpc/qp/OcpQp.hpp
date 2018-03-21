#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/qp/OcpQpBase.hpp>

#include <tmpc/Matrix.hpp>

#include <vector>
#include <initializer_list>

namespace tmpc 
{
	///
	/// \brief One stage of OCP QP.
	///
	/// Solver-agnostic memory layout.
	///
	template <typename Kernel_>
	class OcpQp
	:	public OcpQpBase<OcpQp<Kernel_>>
	{
	public:
		using Kernel = Kernel_;
		using size_type = typename Kernel::size_t;
		using Real = typename Kernel::Real;
		using Matrix = DynamicMatrix<Kernel>;
		using Vector = DynamicVector<Kernel>;
		//using OcpQpBase<OcpQp<Kernel>>::operator=;


		/// @brief Default constructor makes OcpQp DefaultConstructible.
		///
		/// "default" default constructor is fine.
		OcpQp() = default;


		/// @brief Copy constructor makes OcpQp CopyConstructible.
		///
		/// "default" copy constructor is fine.
		OcpQp(OcpQp const&) = default;


		/// @brief Copy assignment makes OcpQp CopyAssignable.
		///
		/// "default" copy assignment is fine.
		/// OcpQp& operator=(OcpQp const&) = default;


		/// @brief Move assignment makes OcpQp MoveAssignable.
		///
		/// "default" move assignment is fine.
		OcpQp& operator=(OcpQp &&) = default;


		OcpQp(OcpSize const& sz, size_type nx_next = 0)
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
		,	Zl_(sz.ns(), sz.ns())
		,	Zu_(sz.ns(), sz.ns())
		,	zl_(sz.ns())
		,	zu_(sz.ns())
		,	idxs_(sz.nx(), 0)
		{
		}


		template <typename Expr>
		OcpQp(OcpQpExpressionBase<Expr> const& expr)
		:	OcpQp(expr.size())
		{
			expr.evalTo(*this);
		}


		template <typename Other>
		OcpQp(OcpQpBase<Other> const& rhs)
		:	size_(rhs.size())
		,	Q_(rhs.Q())
		,	R_(rhs.R())
		,	S_(rhs.S())
		,	q_(rhs.q())
		,	r_(rhs.r())
		,	A_(rhs.A())
		,	B_(rhs.B())
		,	b_(rhs.b())
		,	C_(rhs.C())
		,	D_(rhs.D())
		,	lbx_(rhs.lbx())
		,	ubx_(rhs.ubx())
		,	lbu_(rhs.lbu())
		,	ubu_(rhs.ubu())
		,	lbd_(rhs.lbd())
		,	ubd_(rhs.ubd())
		,	Zl_(rhs.Zl())
		,	Zu_(rhs.Zu())
		,	zl_(rhs.zl())
		,	zu_(rhs.zu())
		,	idxs_(size_.ns())
		{
			this->idxs(rhs.idxs());
		}

		const Matrix& A() const {
			return A_;
		}

		template <typename T>
		void A(const T& a) {
			noresize(A_) = a;
		}

		const Matrix& B() const {
			return B_;
		}

		template <typename T>
		void B(const T& b) {
			noresize(B_) = b;
		}

		Vector const& b() const {
			return b_;
		}

		template <typename T>
		void b(const T& b) {
			noresize(b_) = b;
		}

		const Matrix& C() const {
			return C_;
		}

		template <typename T>
		void C(const T& c) {
			noresize(C_) = c;
		}

		const Matrix& D() const {
			return D_;
		}

		template <typename T>
		void D(const T& d) {
			noresize(D_) = d;
		}

		const Vector& lbd() const {
			return lbd_;
		}

		template <typename T>
		void lbd(const T& lbd) {
			noresize(lbd_) = lbd;
		}

		const Vector& lbu() const {
			return lbu_;
		}

		template <typename T>
		void lbu(const T& lbu) {
			noresize(lbu_) = lbu;
		}

		const Vector& lbx() const {
			return lbx_;
		}

		template <typename T>
		void lbx(const T& lbx) {
			noresize(lbx_) = lbx;
		}

		const Matrix& Q() const {
			return Q_;
		}

		template <typename T>
		void Q(const T& q) {
			noresize(Q_) = q;
		}

		template <typename T>
		auto& Q(std::initializer_list<std::initializer_list<T>> q) 
		{
			noresize(Q_) = q;
			return *this;
		}

		const Matrix& R() const {
			return R_;
		}

		template <typename T>
		void R(const T& r) {
			noresize(R_) = r;
		}

		template <typename T>
		auto& R(std::initializer_list<std::initializer_list<T>> val) 
		{
			noresize(R_) = val;
			return *this;
		}

		const Matrix& S() const {
			return S_;
		}

		template <typename T>
		void S(const T& s) {
			noresize(S_) = s;
		}

		template <typename T>
		auto& S(std::initializer_list<std::initializer_list<T>> val) 
		{
			noresize(S_) = val;
			return *this;
		}

		const Vector& q() const {
			return q_;
		}

		template <typename T>
		void q(const T& q) {
			noresize(q_) = q;
		}

		template <typename T>
		auto& q(std::initializer_list<T> val) 
		{
			noresize(q_) = val;
			return *this;
		}

		const Vector& r() const {
			return r_;
		}

		template <typename T>
		void r(const T& r) {
			noresize(r_) = r;
		}

		template <typename T>
		auto& r(std::initializer_list<T> val) 
		{
			noresize(r_) = val;
			return *this;
		}

		const Vector& ubd() const {
			return ubd_;
		}

		template <typename T>
		void ubd(const T& ubd) {
			noresize(ubd_) = ubd;
		}

		const Vector& ubu() const {
			return ubu_;
		}

		template <typename T>
		void ubu(const T& ubu) {
			noresize(ubu_) = ubu;
		}

		const Vector& ubx() const {
			return ubx_;
		}

		template <typename T>
		void ubx(const T& ubx) {
			noresize(ubx_) = ubx;
		}


		// -----------------------------------------------------------
		// Soft constraints cost
		// -----------------------------------------------------------
		auto const& impl_Zl() const
		{
			return Zl_;
		}


		template <typename T>
		void impl_Zl(T const& val)
		{
			noresize(Zl_) = val;
		}


		auto const& impl_Zu() const
		{
			return Zu_;
		}


		template <typename T>
		void impl_Zu(T const& val)
		{
			noresize(Zu_) = val;
		}


		auto const& impl_zl() const
		{
			return zl_;
		}


		template <typename T>
		void impl_zl(T const& val)
		{
			noresize(zl_) = val;
		}


		auto const& impl_zu() const
		{
			return zu_;
		}


		template <typename T>
		void impl_zu(T const& val)
		{
			noresize(zu_) = val;
		}


		// -----------------------------------------------------------
		// Soft constraints index
		// -----------------------------------------------------------
		auto const& impl_idxs() const
		{
			return idxs_;
		}


		template <typename T> 
		void impl_idxs(T const& val) 
		{ 
			std::copy(val.begin(), val.end(), idxs_.begin());
		}


		// -----------------------------------------------------------
		OcpSize const& size() const
		{
			return size_;
		}

		size_t nxNext() const
		{
			return rows(A_);
		}

	private:
		OcpSize size_;
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
		Matrix Zl_;
		Matrix Zu_;
		Vector zl_;
		Vector zu_;
		std::vector<size_t> idxs_;
	};
}
