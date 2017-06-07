#pragma once

#include "UnsolvedQpException.hpp"
#include <tmpc/Matrix.hpp>
#include <tmpc/qp/QpSize.hpp>

#include <qpOASES.hpp>

#include <ostream>

namespace tmpc {

class QpOasesSolveException : public UnsolvedQpException
{
public:
	QpOasesSolveException(qpOASES::returnValue code);

	qpOASES::returnValue code() const	{ return _code;	}
	char const * what() const noexcept override { return msg_.c_str(); }

private:
	qpOASES::returnValue const _code;
	std::string const msg_;
};

namespace detail {
qpOASES::Options qpOASES_DefaultOptions();
}

/**
 * \brief QP workspace using qpOASES
 */
class QpOasesWorkspace
{
public:
	// Scalar data type.
	using Scalar = double;

	// Matrix storage option -- important!
	// Must be rowMajor, because qpOASES expects input matrices in row-major format.
	static auto constexpr storageOrder = rowMajor;

	typedef unsigned int size_type;
	typedef DynamicMatrix<Scalar, storageOrder> Matrix;
	typedef DynamicVector<Scalar, columnVector> Vector;
	typedef Submatrix<Matrix> SubM;
	typedef Subvector<Vector> SubV;
	typedef CustomVector<Scalar, unaligned, unpadded> VectorMap;

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

		VectorMap const& get_x() const
		{
			return x_;
		}

		VectorMap const& get_u() const
		{
			return u_;
		}

		VectorMap const& get_pi() const
		{
			return pi_;
		}

		VectorMap const& get_lam_u() const
		{
			return lamU_;
		}

		VectorMap const& get_lam_x() const
		{
			return lamX_;
		}

		VectorMap const& get_lam_d() const
		{
			return lam_;
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
		VectorMap x_;
		VectorMap u_;
		VectorMap lamX_;
		VectorMap lamU_;
		VectorMap lam_;
		VectorMap pi_;
	};

	template <typename InputIt>
	QpOasesWorkspace(InputIt sz_begin, InputIt sz_end, qpOASES::Options const& options = detail::qpOASES_DefaultOptions())
	:	QpOasesWorkspace(sz_begin, sz_end, options,
			numVariables(sz_begin, sz_end), numEqualities(sz_begin, sz_end) + numInequalities(sz_begin, sz_end))
	{
		_problem.setOptions(options);
	}

	/**
	 * \brief Copy constructor.
	 *
	 * Copying is not allowed.
	 */
	QpOasesWorkspace(QpOasesWorkspace const&) = delete;

	/**
	 * \brief Move constructor.
	 *
	 * Move-construction is ok.
	 */
	QpOasesWorkspace(QpOasesWorkspace&& rhs) = default;

	QpOasesWorkspace& operator=(QpOasesWorkspace const&) = delete;
	QpOasesWorkspace& operator=(QpOasesWorkspace&&) = delete;

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

	// qpOASES-specific part
	//

	//
	// Full matrix and vector access functions.
	//
	Matrix& H() { return H_; }
	const Matrix& H() const { return H_; }

	Vector& g() { return g_; }
	const Vector& g() const { return g_; }

	Matrix& A() { return A_; }
	const Matrix& A() const { return A_; }

	Vector& lbA() { return lbA_; }
	const Vector& lbA() const { return lbA_; }

	Vector& ubA() { return ubA_; }
	const Vector& ubA() const { return ubA_; }

	Vector& lb() { return lb_; }
	const Vector& lb() const { return lb_; }

	Vector& ub() { return ub_; }
	const Vector& ub() const { return ub_; }

	bool getHotStart() const noexcept { return _hotStart; }

	// Get maximum number of working set recalculations for qpOASES
	unsigned const getMaxWorkingSetRecalculations() const noexcept { return _maxWorkingSetRecalculations; }

	// Set maximum number of working set recalculations for qpOASES
	void setMaxWorkingSetRecalculations(unsigned val) noexcept { _maxWorkingSetRecalculations = val; }

	void solve();

private:
	bool _hotStart = false;

	// TODO: wrap _problem into a pImpl to
	// a) Reduce dependencies
	// b) Avoid deep-copy of qpOASES::SQProblem object of move-construction of CondensingSolver.
	qpOASES::SQProblem _problem;
	unsigned _maxWorkingSetRecalculations = 1000;

	Matrix H_;
	Vector g_;

	// The layout of lb_ is [lbx, lbu, ...]
	Vector lb_;

	// The layout of ub_ is [ubx, ubu, ...]
	Vector ub_;

	Matrix A_;

	// The layout of lbA_ is [lbb, lbd, ...]
	Vector lbA_;

	// The layout of ubA_ is [ubb, ubd, ...]
	Vector ubA_;

	DynamicVector<Scalar> primalSolution_;
	DynamicVector<Scalar> dualSolution_;

	/// \brief Number of iterations performed by the QP solver.
	unsigned numIter_ = 0;

	std::vector<Stage> stage_;

	static auto constexpr sNaN = std::numeric_limits<Scalar>::signaling_NaN();

	template <typename InputIt>
	QpOasesWorkspace(InputIt sz_begin, InputIt sz_end, qpOASES::Options const& options, size_t nx, size_t nc)
	:	_problem(nx, nc)
	,	H_(nx, nx, Scalar{0})
	,	g_(nx, sNaN)
	, 	lb_(nx, sNaN)
	, 	ub_(nx, sNaN)
	, 	A_(nc, nx, Scalar{0})
	, 	lbA_(nc, sNaN)
	, 	ubA_(nc, sNaN)
	,	primalSolution_(nx)
	,	dualSolution_(nx + nc)
	{
		stage_.reserve(std::distance(sz_begin, sz_end));

		double * x = primalSolution_.data();
		double * lambda_x = dualSolution_.data();
		double * lambda = dualSolution_.data() + primalSolution_.size();

		for (auto s = sz_begin; s != sz_end; ++s)
		{
			auto const s_next = s + 1;
			auto const nx_next = s_next != sz_end ? s_next->nx() : 0;
			stage_.emplace_back(*s, nx_next, x, lambda_x, lambda);

			x += s->nx() + s->nu();
			lambda_x += s->nx() + s->nu();
			lambda += s->nc() + nx_next;
		}
		
		_problem.setOptions(options);
	}

	void InitStages();
};

}	// namespace tmpc
