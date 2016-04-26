#include <CondensingSolver.hpp>
#include <qpOASESException.hpp>
#include <Condensing.hpp>

namespace camels
{
	void CondensingSolver::Solve(MultiStageQP const& msqp, Point& solution)
	{
		// Check argument sizes.
		if (!(msqp.nX() == nX() && msqp.nU() == nU() && msqp.nT() == nT()))
			throw std::invalid_argument("CondensingSolver::Solve(): size of MultistageQP does not match solver sizes, sorry.");

		if (!(solution.nX() == nX() && solution.nU() == nU() && solution.nT() == nT()))
			throw std::invalid_argument("CondensingSolver::Solve(): size of solution Point does not match solver sizes, sorry.");

		// Make a condensed problem.
		Condense(msqp, _condensedQP);

		/* Solve the condensed QP. */
		int nWSR = 1000;
		const auto res = _hotStart ?
			_problem.hotstart(_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR) :
			_problem.init    (_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR);

		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw CondensingSolverSolveException(res, _condensedQP);

		_hotStart = true;

		/* Get solution of the condensed QP. */
		_problem.getPrimalSolution(_condensedSolution.data());
		//problem.getDualSolution(yOpt);

		// Calculate the solution of the multi-stage QP.
		solution.w(0).topRows(nX()) = _condensedSolution.topRows(nX());
		for (size_type i = 0; i < nT(); ++i)
		{
			auto z_i = solution.w(i);
			auto x_i = z_i.topRows(nX());
			auto u_i = z_i.bottomRows(nU());
			auto x_i_plus = solution.w(i + 1).topRows(nX());

			u_i = _condensedSolution.middleRows(nX() + i * nU(), nU());
			x_i_plus = msqp.C(i) * z_i + msqp.c(i);
		}
	}


	CondensingSolverSolveException::CondensingSolverSolveException(qpOASES::returnValue code, qpOASESProgram const& cqp) :
		std::runtime_error("CondensingSolver::Solve() failed. qpOASES return code " + std::to_string(code)),
		_code(code), _CondensedQP(cqp)
	{

	}

	qpOASES::returnValue CondensingSolverSolveException::getCode() const
	{
		return _code;
	}

	qpOASESProgram const& CondensingSolverSolveException::getCondensedQP() const
	{
		return _CondensedQP;
	}

	CondensingSolver::Point::Point(size_type nx, size_type nu, size_type nt)
	:	_data((nx + nu) * nt + nx)
	,	_nx(nx)
	,	_nz(nx + nu)
	,	_nt(nt)
	{
	}

	CondensingSolver::Point::VectorMap CondensingSolver::Point::w(unsigned i)
	{
		if(!(i < _nt + 1))
			throw std::out_of_range("CondensingSolver::Point::w(): index is out of range");

		return VectorMap(_data.data() + i * _nz, i < _nt ? _nz : _nx);
	}

	CondensingSolver::Point::VectorConstMap CondensingSolver::Point::w(unsigned i) const
	{
		if (!(i < _nt + 1))
			throw std::out_of_range("CondensingSolver::Point::w(): index is out of range");

		return VectorConstMap(_data.data() + i * _nz, i < _nt ? _nz : _nx);
	}

	void CondensingSolver::Point::shift()
	{
		std::copy_n(_data.begin() + _nz, (_nt - 1) * _nz + _nx, _data.begin());
	}

	CondensingSolver::Point& CondensingSolver::Point::operator+=(Point const& rhs)
	{
		if (rhs.nT() != nT())
			throw std::invalid_argument("CondensingSolver::Point::operator+=(): arguments have different sizes!");

		std::transform(_data.cbegin(), _data.cend(), rhs._data.cbegin(), _data.begin(), std::plus<double>());

		return *this;
	}
}

std::ostream& operator<<(std::ostream& os, camels::CondensingSolver::Point const& point)
{
	for (camels::CondensingSolver::size_type i = 0; i <= point.nT(); ++i)
		os << point.w(i) << std::endl;

	return os;
}
