#pragma once

#include "QpSize.hpp"

#include <vector>
#include <initializer_list>

namespace tmpc 
{
	template <typename Kernel>
	class QuadraticProblemStage
	{
	public:
		using size_type = typename Kernel::size_t;
		using Real = typename Kernel::Real;
		using Matrix = typename Kernel::DynamicMatrix;
		using Vector = typename Kernel::DynamicVector;

		QuadraticProblemStage(QpSize const& sz, size_type nx_next)
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

		template <typename Expr>
		QuadraticProblemStage(Expr const& rhs)
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
		{
		}

		QuadraticProblemStage(QuadraticProblemStage const&) = default; //delete;
		QuadraticProblemStage(QuadraticProblemStage &&) = default;

		template <typename Expr>
		QuadraticProblemStage& operator=(Expr const& rhs)
		{
			assign(*this, rhs);
			return *this;
		}

		const Matrix& A() const {
			return A_;
		}

		template <typename T>
		void A(const T& a) {
			full(A_) = a;
		}

		const Matrix& B() const {
			return B_;
		}

		template <typename T>
		void B(const T& b) {
			full(B_) = b;
		}

		Vector const& b() const {
			return b_;
		}

		template <typename T>
		void b(const T& b) {
			full(b_) = b;
		}

		const Matrix& C() const {
			return C_;
		}

		template <typename T>
		void C(const T& c) {
			full(C_) = c;
		}

		const Matrix& D() const {
			return D_;
		}

		template <typename T>
		void D(const T& d) {
			full(D_) = d;
		}

		const Vector& lbd() const {
			return lbd_;
		}

		template <typename T>
		void lbd(const T& lbd) {
			full(lbd_) = lbd;
		}

		const Vector& lbu() const {
			return lbu_;
		}

		template <typename T>
		void lbu(const T& lbu) {
			full(lbu_) = lbu;
		}

		const Vector& lbx() const {
			return lbx_;
		}

		template <typename T>
		void lbx(const T& lbx) {
			full(lbx_) = lbx;
		}

		const Matrix& Q() const {
			return Q_;
		}

		template <typename T>
		void Q(const T& q) {
			full(Q_) = q;
		}

		const Matrix& R() const {
			return R_;
		}

		template <typename T>
		void R(const T& r) {
			full(R_) = r;
		}

		const Matrix& S() const {
			return S_;
		}

		template <typename T>
		void S(const T& s) {
			full(S_) = s;
		}

		const Vector& q() const {
			return q_;
		}

		template <typename T>
		void q(const T& q) {
			full(q_) = q;
		}

		const Vector& r() const {
			return r_;
		}

		template <typename T>
		void r(const T& r) {
			full(r_) = r;
		}

		const Vector& ubd() const {
			return ubd_;
		}

		template <typename T>
		void ubd(const T& ubd) {
			full(ubd_) = ubd;
		}

		const Vector& ubu() const {
			return ubu_;
		}

		template <typename T>
		void ubu(const T& ubu) {
			full(ubu_) = ubu;
		}

		const Vector& ubx() const {
			return ubx_;
		}

		template <typename T>
		void ubx(const T& ubx) {
			full(ubx_) = ubx;
		}

		QpSize const& size() const
		{
			return size_;
		}

		size_t nxNext() const
		{
			return rows(A_);
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
	template <typename Kernel>
	class QuadraticProblem
	{
	public:
		using size_type = typename Kernel::size_t;
		using Real = typename Kernel::Real;
		using Matrix = typename Kernel::DynamicMatrix;
		using Vector = typename Kernel::DynamicVector;
		using Stage = QuadraticProblemStage<Kernel>;

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
}
