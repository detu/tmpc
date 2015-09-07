#include <CondensingSolver.hpp>
#include <qpOASESException.hpp>

namespace camels
{
	void CondensingSolver::Condense(const MultiStageQP& msqp)
	{
		if (msqp.size() != size())
			throw std::invalid_argument("CondensingSolver::Condense(): the problem has a size different from what the solver expects.");

		Eigen::MatrixXd M(nX(), nIndep());
		M.setZero();

		Eigen::VectorXd v(nX());
		v.setZero();

		auto& Hc = _condensedQP.H();
		auto& gc = _condensedQP.g();

		Hc.setZero();
		gc.setZero();

		_condensedQP.lb().topRows(nX()) = msqp.xMin(0);
		_condensedQP.ub().topRows(nX()) = msqp.xMax(0);

		for (unsigned k = 0; k <= nT(); ++k)
		{
			auto M_k = M.leftCols(nX() + k * nU());
			if (k == 0)
			{
				M_k.setIdentity();
				v.setZero();
			}
			else
			{
				const auto A_k_minus = msqp.C(k - 1).leftCols(nX());
				const auto B_k_minus = msqp.C(k - 1).rightCols(nU());
				M_k = A_k_minus * M_k;
				M_k.rightCols(nU()) = B_k_minus;
				v = A_k_minus * v + msqp.c(k - 1);

				_condensedQP.A()  .middleRows((k - 1) * nX(), nX()) = M;
				_condensedQP.lbA().middleRows((k - 1) * nX(), nX()) = msqp.xMin(k) - v;
				_condensedQP.ubA().middleRows((k - 1) * nX(), nX()) = msqp.xMax(k) - v;
			}

			const auto H_k = msqp.H(k);
			const auto g_k = msqp.g(k);

			if (k < nT())
			{
				const auto Q = H_k.topLeftCorner(nX(), nX()).selfadjointView<Eigen::Upper>();
				const auto S = H_k.topRightCorner(nX(), nU());
				const auto R = H_k.bottomRightCorner(nU(), nU());

				const auto nn = M_k.cols();
				auto Hc_k = Hc.topLeftCorner(nn + nU(), nn + nU());
				Hc_k.topLeftCorner(nn, nn).triangularView<Eigen::Upper>() += M_k.transpose() * Q * M_k;
				Hc_k.topRightCorner(nn, nU()) += M_k.transpose() * S;
				Hc_k.bottomRightCorner(nU(), nU()).triangularView<Eigen::Upper>() += R;

				auto gc_k = gc.topRows(nn + nU());
				gc_k.topRows(nn) += M_k.transpose() * (g_k.topRows(nX()) + Q * v);
				gc_k.bottomRows(nU()) += g_k.bottomRows(nU()) + S.transpose() * v;

				_condensedQP.lb().middleRows(nX() + k * nU(), nU()) = msqp.uMin(k);
				_condensedQP.ub().middleRows(nX() + k * nU(), nU()) = msqp.uMax(k);
			}
			else
			{
				// Final state.
				Hc.triangularView<Eigen::Upper>() += M_k.transpose() * H_k.selfadjointView<Eigen::Upper>() * M_k;
				gc += M_k.transpose() * (g_k.topRows(nX()) + H_k * v);
			}
		}

		Hc = Hc.selfadjointView<Eigen::Upper>();
	}

	void CondensingSolver::Solve(const MultiStageQP& msqp)
	{
		// Make a condensed problem.
		Condense(msqp);

		/* Solve the condensed QP. */
		int nWSR = 1000;
		const auto res = _hotStart ?
			_problem.hotstart(_condensedQP.H().data(), _condensedQP.g().data(), _condensedQP.A().data(),
			_condensedQP.lb().data(), _condensedQP.ub().data(), _condensedQP.lbA().data(), _condensedQP.ubA().data(), nWSR) :
			_problem.init(_condensedQP.H().data(), _condensedQP.g().data(), _condensedQP.A().data(),
			_condensedQP.lb().data(), _condensedQP.ub().data(), _condensedQP.lbA().data(), _condensedQP.ubA().data(), nWSR);

		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw qpOASES::Exception(res);

		_hotStart = true;

		/* Get solution of the condensed QP. */
		_problem.getPrimalSolution(_primalCondensedSolution.data());
		//problem.getDualSolution(yOpt);

		// Calculate the solution of the multi-stage QP.
		_primalSolution.topRows(nX()) = _primalCondensedSolution.topRows(nX());
		for (size_type i = 0; i < nT(); ++i)
		{
			auto z_i = _primalSolution.middleRows(i * nZ(), nZ());
			auto x_i = z_i.topRows(nX());
			auto u_i = z_i.bottomRows(nU());
			auto x_i_plus = _primalSolution.middleRows((i + 1) * nZ(), nX());
			
			u_i = _primalCondensedSolution.middleRows(nX() + i * nU(), nU());
			x_i_plus = msqp.C(i) * z_i + msqp.c(i);
		}
	}

	const CondensingSolver::Vector& CondensingSolver::getPrimalSolution() const
	{
		return _primalSolution;
	}

	CondensingSolver::CondensingSolver(size_type nx, size_type nu, size_type nt) :
		CondensingSolver(MultiStageQPSize(nx, nu, 0, 0, nt))
	{		
	}

	CondensingSolver::CondensingSolver(const MultiStageQPSize& size) :
		_condensedQP(size.nIndep(), size.nDep()),
		_size(size),
		_primalCondensedSolution(size.nIndep()),
		_primalSolution(size.nVar()),
		_problem(size.nIndep(), size.nDep())
	{
		qpOASES::Options options;
		_problem.setOptions(options);
	}

	const CondensingSolver::Vector& CondensingSolver::getPrimalCondensedSolution() const
	{
		return _primalCondensedSolution;
	}

	camels::CondensingSolver::size_type CondensingSolver::nT() const
	{
		return _size.nT();
	}

	camels::CondensingSolver::size_type CondensingSolver::nX() const
	{
		return _size.nX();
	}

	camels::CondensingSolver::size_type CondensingSolver::nZ() const
	{
		return _size.nZ();
	}

	camels::CondensingSolver::size_type CondensingSolver::nU() const
	{
		return _size.nU();
	}

	camels::CondensingSolver::size_type CondensingSolver::nD() const
	{
		return _size.nD();
	}

	camels::CondensingSolver::size_type CondensingSolver::nDT() const
	{
		return _size.nDT();
	}

	camels::CondensingSolver::size_type CondensingSolver::nIndep() const
	{
		return _size.nIndep();
	}

	camels::CondensingSolver::size_type CondensingSolver::nDep() const
	{
		return _size.nDep();
	}

	camels::CondensingSolver::size_type CondensingSolver::nVar() const
	{
		return _size.nVar();
	}

	const MultiStageQPSize& CondensingSolver::size() const
	{
		return _size;
	}
}
