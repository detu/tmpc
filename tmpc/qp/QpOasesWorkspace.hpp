#pragma once

#include "QpSolverException.hpp"

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>

#include <qpOASES.hpp>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <ostream>
#include <vector>
#include <memory>

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

private:

	using size_t = typename Kernel::size_t;

	class Stage
	:	public OcpSolutionBase<Stage>
	,	public OcpQpBase<Stage>
	{
	public:
		using Kernel = QpOasesWorkspace::Kernel;

		Stage(Workspace& ws, OcpSize const& sz, size_t n, size_t na, size_t nx_next)
		:	ws_(&ws)
		,	n_(n)
		,	na_(na)
		,	size_(sz)
		,	nxNext_(nx_next)
		{
			this->setNaN();
		}

		Stage(Stage const&) = delete;
		Stage(Stage && s) = default;

		// Set new workspace.
		void workspace(Workspace& ws)
		{
			ws_ = &ws;
		}

		auto Q() const 
		{	
			return submatrix(ws_->H, n_, n_, size_.nx(), size_.nx());
		}
		
		template <typename T> 
		void Q(const T& q) 
		{ 
			submatrix(ws_->H, n_, n_, size_.nx(), size_.nx()) = q; 
		}

		auto R() const 
		{	
			return submatrix(ws_->H, n_ + size_.nx(), n_ + size_.nx(), size_.nu(), size_.nu()); 
		}

		template <typename T> 
		void R(const T& r) 
		{ 
			submatrix(ws_->H, n_ + size_.nx(), n_ + size_.nx(), size_.nu(), size_.nu()) = r; 
		}

		auto S() const 
		{	
			return submatrix(ws_->H, n_, n_ + size_.nx(), size_.nx(), size_.nu()); 
		}

		template <typename T> 
		void S(const T& s) 
		{
			// Set both S and S^T in H
			submatrix(ws_->H, n_, n_ + size_.nx(), size_.nx(), size_.nu()) = s;
			submatrix(ws_->H, n_ + size_.nx(), n_, size_.nu(), size_.nx()) = trans(s); 
		}

		void S(Real val) 
		{
			// Set both S and S^T in H
			submatrix(ws_->H, n_, n_ + size_.nx(), size_.nx(), size_.nu()) = val;
			submatrix(ws_->H, n_ + size_.nx(), n_, size_.nu(), size_.nx()) = val; 
		}

		auto q() const 
		{	
			return subvector(ws_->g, n_, size_.nx());
		}

		template <typename T> 
		void q(const T& q) 
		{ 
			subvector(ws_->g, n_, size_.nx()) = q; 
		}

		auto r() const 
		{	
			return subvector(ws_->g , n_ + size_.nx(), size_.nu());
		}

		template <typename T> 
		void r(const T& r) 
		{ 
			subvector(ws_->g , n_ + size_.nx(), size_.nu()) = r;
		}

		auto A() const 
		{	
			return submatrix(ws_->A, na_, n_, nxNext_, size_.nx()); 
		}

		template <typename T> 
		void A(const T& a) 
		{ 
			submatrix(ws_->A, na_, n_, nxNext_, size_.nx()) = a; 
		}

		auto B() const 
		{	
			return submatrix(ws_->A, na_, n_ + size_.nx(), nxNext_, size_.nu()); 
		}

		template <typename T> 
		void B(const T& b) 
		{ 
			submatrix(ws_->A, na_, n_ + size_.nx(), nxNext_, size_.nu()) = b; 
		}

		auto b() const 
		{ 
			return -subvector(ws_->lbA, na_, nxNext_); 
		}

		template <typename T> 
		void b(const T& b) 
		{
			// Set lbA = ubA = -b
			subvector(ws_->lbA, na_, nxNext_) = subvector(ws_->ubA, na_, nxNext_) = -b; 
		}

		auto C() const 
		{	
			return submatrix(ws_->A, na_ + nxNext_, n_, size_.nc(), size_.nx()); 
		}

		template <typename T> 
		void C(const T& c) 
		{ 
			submatrix(ws_->A, na_ + nxNext_, n_, size_.nc(), size_.nx()) = c; 
		}

		auto D() const 
		{	
			return submatrix(ws_->A, na_ + nxNext_, n_ + size_.nx(), size_.nc(), size_.nu()); 
		}

		template <typename T> 
		void D(const T& d) 
		{ 
			submatrix(ws_->A, na_ + nxNext_, n_ + size_.nx(), size_.nc(), size_.nu()) = d;
		}

		auto lbd() const 
		{ 
			return subvector(ws_->lbA, na_ + nxNext_, size_.nc()); 
		}

		template <typename T> 
		void lbd(const T& lbd) 
		{ 
			subvector(ws_->lbA, na_ + nxNext_, size_.nc()) = lbd;
		}

		auto ubd() const 
		{ 
			return subvector(ws_->ubA, na_ + nxNext_, size_.nc()); 
		}

		template <typename T> 
		void ubd(const T& ubd) 
		{ 
			subvector(ws_->ubA, na_ + nxNext_, size_.nc()) = ubd; 
		}

		auto lbu() const 
		{ 
			return subvector(ws_->lb, n_ + size_.nx(), size_.nu()); 
		}

		template <typename T> 
		void lbu(const T& lbu) 
		{ 
			subvector(ws_->lb, n_ + size_.nx(), size_.nu()) = lbu; 
		}

		auto ubu() const 
		{ 
			return subvector(ws_->ub, n_ + size_.nx(), size_.nu()); 
		}

		template <typename T> 
		void ubu(const T& ubu) 
		{ 
			subvector(ws_->ub, n_ + size_.nx(), size_.nu()) = ubu; 
		}

		auto lbx() const 
		{ 
			return subvector(ws_->lb, n_, size_.nx()); 
		}

		template <typename T> 
		void lbx(const T& lbx) 
		{ 
			subvector(ws_->lb, n_, size_.nx()) = lbx; 
		}		

		auto ubx() const 
		{ 
			return subvector(ws_->ub, n_, size_.nx()); 
		}

		template <typename T> 
		void ubx(const T& ubx) 
		{ 
			subvector(ws_->ub, n_, size_.nx()) = ubx; 
		}

		auto x() const 
		{ 
			return subvector(ws_->primalSolution, n_, size_.nx());	
		}

		auto u() const 
		{ 
			return subvector(ws_->primalSolution, n_ + size_.nx(), size_.nu());	
		}

		auto pi() const	
		{ 
			return subvector(ws_->dualSolution, ws_->primalSolution.size() + na_ + size_.nc(), nxNext_);
		}

		auto lam_lbu() const 
		{ 
			// lamU_(subvector(ws_->dualSolution, n + size_.nx(), size_.nu()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NU.
			return subvector(ws_->dualSolution, n_ + size_.nx(), size_.nu());
		}

		auto lam_ubu() const 
		{ 
			// lamU_(subvector(ws_->dualSolution, n + size_.nx(), size_.nu()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NU.
			return subvector(ws_->dualSolution, n_ + size_.nx(), size_.nu()); 
		}

		auto lam_lbx() const 
		{ 
			// lamX_(subvector(ws_->dualSolution, n          , size_.nx()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NX.
			return subvector(ws_->dualSolution, n_, size_.nx()); 
		}

		auto lam_ubx() const 
		{ 
			// lamX_(subvector(ws_->dualSolution, n          , size_.nx()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NX.
			return subvector(ws_->dualSolution, n_, size_.nx()); 
		}

		auto lam_lbd() const 
		{ 
			// lam_(subvector(ws_->dualSolution, ws_->primalSolution.size() + na, size_.nc()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NC.
			return subvector(ws_->dualSolution, ws_->primalSolution.size() + na_, size_.nc()); 
		}

		auto lam_ubd() const 
		{ 
			// lam_(subvector(ws_->dualSolution, ws_->primalSolution.size() + na, size_.nc()))

			// TODO: this must be fixed, since qpOASES returns both lower an upper bound 
			// Lagrange multipliers in one vector of size NC.
			return subvector(ws_->dualSolution, ws_->primalSolution.size() + na_, size_.nc()); 
		}

		OcpSize const& size() const { return size_; }

	private:
		Workspace * ws_;
		size_t const n_;
		size_t const na_;
		OcpSize const size_;
		size_t const nxNext_;
	};

public:
	template <typename InputIt>
	QpOasesWorkspace(InputIt sz_begin, InputIt sz_end)
	:	problem_(numVariables(sz_begin, sz_end), numEqualities(sz_begin, sz_end) + numInequalities(sz_begin, sz_end))
	,	ws_(problem_.getNV(), problem_.getNC())
	{
		ws_.initABlocks(sz_begin, sz_end);
		initStages(sz_begin, sz_end);
		
		problem_.setOptions(detail::qpOASES_DefaultOptions());
	}

	explicit QpOasesWorkspace(std::initializer_list<OcpSize> sz)
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
	 */
	QpOasesWorkspace(QpOasesWorkspace && rhs)
	:	_hotStart(rhs._hotStart)
	,	problem_(std::move(rhs.problem_))
	,	_maxWorkingSetRecalculations(rhs._maxWorkingSetRecalculations)
	,	ws_(std::move(rhs.ws_))
	,	numIter_(rhs.numIter_)
	,	stage_(std::move(rhs.stage_))
	{
		// Set workspace pointers for the stages.
		for (auto& s : stage_)
			s.workspace(ws_);
	}

	QpOasesWorkspace& operator=(QpOasesWorkspace const&) = delete;
	QpOasesWorkspace& operator=(QpOasesWorkspace &&) = delete;

	class ProblemIterator
	:	public boost::iterator_adaptor<
			ProblemIterator	// derived
		,	typename std::vector<Stage>::iterator	// base
		,	OcpQpBase<Stage>&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	OcpQpBase<Stage>&	// reference
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
		,	OcpQpBase<Stage> const&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	OcpQpBase<Stage> const&	// reference
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
		,	OcpSolutionBase<Stage> const&	// value
		,	boost::random_access_traversal_tag	// category of traversal
		,	OcpSolutionBase<Stage> const&	// reference
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
	auto& H() { return ws_.H; }
	const auto& H() const { return ws_.H; }

	auto& g() { return ws_.g; }
	const auto& g() const { return ws_.g; }

	auto& A() { return ws_.A; }
	const auto& A() const { return ws_.A; }

	auto& lbA() { return ws_.lbA; }
	const auto& lbA() const { return ws_.lbA; }

	auto& ubA() { return ws_.ubA; }
	const auto& ubA() const { return ws_.ubA; }

	auto& lb() { return ws_.lb; }
	const auto& lb() const { return ws_.lb; }

	auto& ub() { return ws_.ub; }
	const auto& ub() const { return ws_.ub; }

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
	// Matrix storage option -- important!
	// Must be rowMajor, because qpOASES expects input matrices in row-major format.
	static auto constexpr storageOrder = rowMajor;

	bool _hotStart = false;

	// TODO (??): wrap problem_ into a pImpl to
	// a) Reduce dependencies
	// b) Avoid deep-copy of qpOASES::SQProblem object on move-construction.
	qpOASES::SQProblem problem_;
	unsigned _maxWorkingSetRecalculations = 1000;
	
	struct Workspace
	{
		Workspace(size_t nx, size_t nc)
		:	memH_(new Real[nx * nx])
		,	H {memH_.get(), nx, nx}
		,	g(nx, sNaN<Real>())
		, 	lb(nx, sNaN<Real>())
		, 	ub(nx, sNaN<Real>())
		,	memA_(new Real[nc * nx])
		, 	A {memA_.get(), nc, nx}
		, 	lbA(nc, Real {})
		, 	ubA(nc, Real {})
		,	primalSolution(nx)
		,	dualSolution(nx + nc)
		{		
			H = Real {0};
			A = Real {0};
		}

		//
		// Init the -I blocks in the A matrix and 0 blocks in lbA and ubA
		//
		template <typename InputIt>
		void initABlocks(InputIt sz_begin, InputIt sz_end)
		{
			// (i, j) = top left corner of the current AB block
			size_t i = 0;
			size_t j = 0;

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
		
		// Memory array for H.
		std::unique_ptr<Real[]> memH_;

		CustomMatrix<Kernel, unaligned, unpadded, storageOrder> H;
		DynamicVector<Kernel> g;

		// The layout of lb is [lbx, lbu, ...]
		DynamicVector<Kernel> lb;

		// The layout of ub is [ubx, ubu, ...]
		DynamicVector<Kernel> ub;

		// Memory array for A.
		std::unique_ptr<Real[]> memA_;

		CustomMatrix<Kernel, unaligned, unpadded, storageOrder> A;

		// The layout of lbA is [lbb, lbd, ...]
		DynamicVector<Kernel> lbA;

		// The layout of ubA_ is [ubb, ubd, ...]
		DynamicVector<Kernel> ubA;

		DynamicVector<Kernel> primalSolution;
		DynamicVector<Kernel> dualSolution;
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
