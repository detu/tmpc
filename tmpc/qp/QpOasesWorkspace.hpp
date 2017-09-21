#pragma once

#include "QpSolverException.hpp"

#include <tmpc/Matrix.hpp>
#include <tmpc/qp/QpSize.hpp>
#include <tmpc/qp/QpStageSolutionBase.hpp>
#include <tmpc/qp/QpStageBase.hpp>

#include <qpOASES.hpp>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <ostream>
#include <vector>

namespace tmpc {

class QpOasesException : public QpSolverException
{
public:
	QpOasesException(qpOASES::returnValue code);

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
template <typename Kernel_>
class QpOasesWorkspace
{
	struct Workspace;

public:
	// Real data type.
	using Real = double;
	
	// Kernel type
	using Kernel = Kernel_;

	// Matrix storage option -- important!
	// Must be rowMajor, because qpOASES expects input matrices in row-major format.
	static auto constexpr storageOrder = rowMajor;

	typedef unsigned int size_type;
	typedef DynamicMatrix<Kernel, storageOrder> Matrix;
	typedef DynamicVector<Kernel, columnVector> Vector;

	class Stage
	:	public QpStageSolutionBase<Stage>
	,	public QpStageBase<Stage>
	{
		typedef Submatrix<Kernel, Matrix> SubM;
		typedef Subvector<Kernel, Vector> SubV;

	public:
		Stage(Workspace& ws, QpSize const& sz, size_t n, size_t na, size_t nx_next)
		:	size_(sz)
		,	Q_ (submatrix(ws.H, n          , n          , sz.nx(), sz.nx()))
		,	R_ (submatrix(ws.H, n + sz.nx(), n + sz.nx(), sz.nu(), sz.nu()))
		,	S_ (submatrix(ws.H, n          , n + sz.nx(), sz.nx(), sz.nu()))
		,	ST_(submatrix(ws.H, n + sz.nx(), n          , sz.nu(), sz.nx()))
		,	q_  (subvector(ws.g , n          , sz.nx()))
		,	r_  (subvector(ws.g , n + sz.nx(), sz.nu()))
		,	lbx_(subvector(ws.lb, n          , sz.nx()))
		,	ubx_(subvector(ws.ub, n          , sz.nx()))
		,	lbu_(subvector(ws.lb, n + sz.nx(), sz.nu()))
		,	ubu_(subvector(ws.ub, n + sz.nx(), sz.nu()))
		,	A_(submatrix(ws.A, na          , n          , nx_next, sz.nx()))
		,	B_(submatrix(ws.A, na          , n + sz.nx(), nx_next, sz.nu()))
		,	lbb_(subvector(ws.lbA, na, nx_next))
		,	ubb_(subvector(ws.ubA, na, nx_next))
		,	C_(submatrix(ws.A, na + nx_next, n          , sz.nc(), sz.nx()))
		,	D_(submatrix(ws.A, na + nx_next, n + sz.nx(), sz.nc(), sz.nu()))
		,	lbd_(subvector(ws.lbA, na + nx_next, sz.nc()))
		,	ubd_(subvector(ws.ubA, na + nx_next, sz.nc()))
		,	x_(subvector(ws.primalSolution, n          , sz.nx()))
		,	u_(subvector(ws.primalSolution, n + sz.nx(), sz.nu()))
		,	lamX_(subvector(ws.dualSolution, n          , sz.nx()))
		,	lamU_(subvector(ws.dualSolution, n + sz.nx(), sz.nu()))
		,	lam_(subvector(ws.dualSolution, ws.primalSolution.size() + na          , sz.nc()))
		,	pi_ (subvector(ws.dualSolution, ws.primalSolution.size() + na + sz.nc(), nx_next))
		{
		}

		Stage(Stage const&) = delete;
		Stage(Stage && s) = default;

		template <typename Expr>
		Stage& operator=(Expr const& rhs)
		{
			assign(*this, rhs);
			return *this;
		}

		const SubM& Q() const {	return Q_; }
		template <typename T> void Q(const T& q) { Q_ = q; }

		const SubM& R() const {	return R_; }
		template <typename T> void R(const T& r) { R_ = r; }

		const SubM& S() const {	return S_; }
		template <typename T> void S(const T& s) { ST_ = (S_ = s).transpose(); }

		const SubV& q() const {	return q_; }
		template <typename T> void q(const T& q) { q_ = q; }

		const SubV& r() const {	return r_; }
		template <typename T> void r(const T& r) { r_ = r; }

		const SubM& A() const {	return A_; }
		template <typename T> void A(const T& a) { A_ = a; }

		const SubM& B() const {	return B_; }
		template <typename T> void B(const T& b) { B_ = b; }

		decltype(auto) b() const { return -lbb_; }
		template <typename T> void b(const T& b) { lbb_ = ubb_ = -b; }

		const SubM& C() const {	return C_; }
		template <typename T> void C(const T& c) { C_ = c; }

		const SubM& D() const {	return D_; }
		template <typename T> void D(const T& d) { D_ = d; }

		const SubV& lbd() const { return lbd_; }
		template <typename T> void lbd(const T& lbd) { lbd_ = lbd; }

		const SubV& ubd() const { return ubd_; }
		template <typename T> void ubd(const T& ubd) { ubd_ = ubd; }

		const SubV& lbu() const { return lbu_; }
		template <typename T> void lbu(const T& lbu) { lbu_ = lbu; }

		const SubV& ubu() const { return ubu_; }
		template <typename T> void ubu(const T& ubu) { ubu_ = ubu; }

		const SubV& lbx() const { return lbx_; }
		template <typename T> void lbx(const T& lbx) { lbx_ = lbx; }		

		const SubV& ubx() const { return ubx_; }
		template <typename T> void ubx(const T& ubx) { ubx_ = ubx; }

		SubV const& x() const { return x_;	}
		SubV const& u() const { return u_;	}
		SubV const& pi() const	{ return pi_; }
		decltype(auto) lam_lbu() const { return subvector(lamU_, 0, size_.nu()); }
		decltype(auto) lam_ubu() const { return subvector(lamU_, size_.nu(), size_.nu()); }
		decltype(auto) lam_lbx() const { return subvector(lamX_, 0, size_.nx()); }
		decltype(auto) lam_ubx() const { return subvector(lamX_, size_.nx(), size_.nx()); }
		decltype(auto) lam_lbd() const { return subvector(lam_, 0, size_.nc()); }
		decltype(auto) lam_ubd() const { return subvector(lam_, size_.nc(), size_.nc()); }

		QpSize const& size() const { return size_; }

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
		SubV x_;
		SubV u_;
		SubV lamX_;
		SubV lamU_;
		SubV lam_;
		SubV pi_;
	};

	template <typename InputIt>
	QpOasesWorkspace(InputIt sz_begin, InputIt sz_end)
	:	problem_(numVariables(sz_begin, sz_end), numEqualities(sz_begin, sz_end) + numInequalities(sz_begin, sz_end))
	,	ws_(problem_.getNV(), problem_.getNC())
	{
		ws_.initABlocks(sz_begin, sz_end);
		initStages(sz_begin, sz_end);
		
		problem_.setOptions(detail::qpOASES_DefaultOptions());
	}

	explicit QpOasesWorkspace(std::initializer_list<QpSize> sz)
	:	QpOasesWorkspace(sz.begin(), sz.end())
	{
	}

	template <typename IteratorRange>
	explicit QpOasesWorkspace(IteratorRange sz)
	:	QpOasesWorkspace(sz.begin(), sz.end())
	{
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
	QpOasesWorkspace(QpOasesWorkspace &&) = default;

	QpOasesWorkspace& operator=(QpOasesWorkspace const&) = delete;
	QpOasesWorkspace& operator=(QpOasesWorkspace &&) = delete;

	class ProblemIterator
	:	public boost::iterator_adaptor<
			ProblemIterator	// derived
		,	typename std::vector<Stage>::iterator	// base
		,	QpStageBase<Stage>&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	QpStageBase<Stage>&	// reference
		>
	{
	public:
		ProblemIterator() = default;
		
		ProblemIterator(typename ProblemIterator::base_type const& p)
		:	ProblemIterator::iterator_adaptor_(p)
		{			
		}
	};

	class ConstProblemIterator
	:	public boost::iterator_adaptor<
			ConstProblemIterator	// derived
		,	typename std::vector<Stage>::const_iterator	// base
		,	QpStageBase<Stage> const&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	QpStageBase<Stage> const&	// reference
		>
	{
	public:
		ConstProblemIterator() = default;
		
		ConstProblemIterator(typename ConstProblemIterator::base_type const& p)
		:	ConstProblemIterator::iterator_adaptor_(p)
		{			
		}
	};

	class ConstSolutionIterator
	:	public boost::iterator_adaptor<
			ConstSolutionIterator	// derived
		,	typename std::vector<Stage>::const_iterator	// base
		,	QpStageSolutionBase<Stage> const&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	QpStageSolutionBase<Stage> const&	// reference
		>
	{
	public:
		ConstSolutionIterator() = default;
		
		ConstSolutionIterator(typename ConstSolutionIterator::base_type const& p)
		:	ConstSolutionIterator::iterator_adaptor_(p)
		{			
		}
	};

	boost::iterator_range<ProblemIterator> problem()
	{
		return boost::iterator_range<ProblemIterator>(stage_.begin(), stage_.end());
	}

	boost::iterator_range<ConstProblemIterator> problem() const
	{
		return boost::iterator_range<ConstProblemIterator>(stage_.begin(), stage_.end());
	}

	boost::iterator_range<ConstSolutionIterator> solution() const
	{
		return boost::iterator_range<ConstSolutionIterator>(stage_.begin(), stage_.end());
	}

	// qpOASES-specific part
	//

	//
	// Full matrix and vector access functions.
	//
	Matrix& H() { return ws_.H; }
	const Matrix& H() const { return ws_.H; }

	Vector& g() { return ws_.g; }
	const Vector& g() const { return ws_.g; }

	Matrix& A() { return ws_.A; }
	const Matrix& A() const { return ws_.A; }

	Vector& lbA() { return ws_.lbA; }
	const Vector& lbA() const { return ws_.lbA; }

	Vector& ubA() { return ws_.ubA; }
	const Vector& ubA() const { return ws_.ubA; }

	Vector& lb() { return ws_.lb; }
	const Vector& lb() const { return ws_.lb; }

	Vector& ub() { return ws_.ub; }
	const Vector& ub() const { return ws_.ub; }

	bool hotStart() const noexcept { return _hotStart; }

	// Get maximum number of working set recalculations for qpOASES
	unsigned const maxWorkingSetRecalculations() const noexcept { return _maxWorkingSetRecalculations; }

	// Set maximum number of working set recalculations for qpOASES
	void maxWorkingSetRecalculations(unsigned val) noexcept { _maxWorkingSetRecalculations = val; }

	/// \brief Number of working set recalculations during last solve().
	unsigned numIter() const { return numIter_; }

	void solve()
	{
		/* Solve the QP. */
		int nWSR = static_cast<int>(_maxWorkingSetRecalculations);
		const auto res = _hotStart ?
			problem_.hotstart(ws_.H.data(), ws_.g.data(), ws_.A.data(),
					ws_.lb.data(), ws_.ub.data(), ws_.lbA.data(), ws_.ubA.data(), nWSR) :
			problem_.init    (ws_.H.data(), ws_.g.data(), ws_.A.data(),
					ws_.lb.data(), ws_.ub.data(), ws_.lbA.data(), ws_.ubA.data(), nWSR);
	
		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw QpOasesException(res);
	
		numIter_ = nWSR;
		_hotStart = true;
	
		/* Get solution data. */
		problem_.getPrimalSolution(ws_.primalSolution.data());
		problem_.getDualSolution(ws_.dualSolution.data());
	}

private:
	bool _hotStart = false;

	// TODO: wrap problem_ into a pImpl to
	// a) Reduce dependencies
	// b) Avoid deep-copy of qpOASES::SQProblem object on move-construction of CondensingSolver.
	qpOASES::SQProblem problem_;
	unsigned _maxWorkingSetRecalculations = 1000;
	
	struct Workspace
	{
		Workspace(size_t nx, size_t nc)
		:	H(nx, nx, Real{0})
		,	g(nx, std::numeric_limits<Real>::signaling_NaN())
		, 	lb(nx, std::numeric_limits<Real>::signaling_NaN())
		, 	ub(nx, std::numeric_limits<Real>::signaling_NaN())
		, 	A(nc, nx, Real{0})
		, 	lbA(nc, std::numeric_limits<Real>::signaling_NaN())
		, 	ubA(nc, std::numeric_limits<Real>::signaling_NaN())
		,	primalSolution(nx)
		,	dualSolution(nx + nc)
		{		
		}

		//
		// Init the -I blocks in the A matrix and 0 blocks in lbA and ubA
		//
		template <typename InputIt>
		void initABlocks(InputIt sz_begin, InputIt sz_end)
		{
			// (i, j) = top left corner of the current AB block
			size_type i = 0;
			size_type j = 0;

			for (auto sz = sz_begin; sz + 1 < sz_end; ++sz)
			{
				// Move (i, j) one column right from the top left corner of the current AB block,
				// which is the top left corner of the -I block.
				j += sz->nx() + sz->nu();

				// Size of the -I block
				auto const nx_next = (sz + 1)->nx();

				// Assign the -I block in A
				submatrix(A, i, j, nx_next, nx_next) = -IdentityMatrix<Kernel>(nx_next);

				// Assign the 0 blocks in lbA and ubA
				subvector(lbA, i, nx_next) = 0.;
				subvector(ubA, i, nx_next) = 0.;

				// Move (i, j) to the top left corner of the next AB block.
				i += nx_next + sz->nc();
			}
		}
		
		Matrix H;
		Vector g;

		// The layout of lb is [lbx, lbu, ...]
		Vector lb;

		// The layout of ub is [ubx, ubu, ...]
		Vector ub;

		Matrix A;

		// The layout of lbA is [lbb, lbd, ...]
		Vector lbA;

		// The layout of ubA_ is [ubb, ubd, ...]
		Vector ubA;

		Vector primalSolution;
		Vector dualSolution;
	};

	Workspace ws_;

	/// \brief Number of iterations performed by the QP solver on the last call to solve().
	unsigned numIter_ = 0;

	std::vector<Stage> stage_;

	//
	// Init stage objects.
	//
	template <typename InputIt>
	void initStages(InputIt sz_begin, InputIt sz_end)
	{
		stage_.reserve(std::distance(sz_begin, sz_end));

		std::size_t n = 0, na = 0;

		for (auto sz = sz_begin; sz != sz_end; ++sz)
		{
			auto const nx_next = sz + 1 != sz_end ? (sz + 1)->nx() : 0;
			stage_.emplace_back(ws_, *sz, n, na, nx_next);

			n += sz->nx() + sz->nu();
			na += nx_next + sz->nc();
		}
	}
};

}	// namespace tmpc
