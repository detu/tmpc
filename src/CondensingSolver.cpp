#include <CondensingSolver.hpp>

namespace camels
{
	struct qpOASESException : std::runtime_error
	{
		qpOASESException(qpOASES::returnValue ret)
		: std::runtime_error(ErrorMessage(ret)), _returnValue(ret)
		{}

		qpOASES::returnValue getReturnValue() const { return _returnValue; }

	private:
		qpOASES::returnValue _returnValue;

		static std::string ErrorMessage(qpOASES::returnValue ret)
		{
			return "qpOASES return code " + std::to_string(ret);
		}
	};

	void CondensingSolver::Condense(const MultiStageQP& msqp)
	{
		assert(msqp.nX() == _Nx && msqp.nU() == _Nu && msqp.nT() == _Nt);

		Eigen::MatrixXd M(_Nx, nIndep());
		M.setZero();

		Eigen::VectorXd v(_Nx);
		v.setZero();

		auto& Hc = _condensedQP.H();
		auto& gc = _condensedQP.g();

		Hc.setZero();
		gc.setZero();

		_condensedQP.lb().topRows(_Nx) = msqp.xMin(0);
		_condensedQP.ub().topRows(_Nx) = msqp.xMax(0);

		for (unsigned k = 0; k <= _Nt; ++k)
		{
			auto M_k = M.leftCols(_Nx + k * _Nu);
			if (k == 0)
			{
				M_k.setIdentity();
				v.setZero();
			}
			else
			{
				const auto A_k_minus = msqp.C(k - 1).leftCols(_Nx);
				const auto B_k_minus = msqp.C(k - 1).rightCols(_Nu);
				M_k = A_k_minus * M_k;
				M_k.rightCols(_Nu) = B_k_minus;
				v = A_k_minus * v + msqp.c(k - 1);

				_condensedQP.A()  .middleRows((k - 1) * _Nx, _Nx) = M;
				_condensedQP.lbA().middleRows((k - 1) * _Nx, _Nx) = msqp.xMin(k) - v;
				_condensedQP.ubA().middleRows((k - 1) * _Nx, _Nx) = msqp.xMax(k) - v;
			}

			const auto H_k = msqp.H(k);
			const auto g_k = msqp.g(k);

			if (k < _Nt)
			{
				const auto Q = H_k.topLeftCorner(_Nx, _Nx);
				const auto S = H_k.topRightCorner(_Nx, _Nu);
				const auto ST = H_k.bottomLeftCorner(_Nu, _Nx);
				const auto R = H_k.bottomRightCorner(_Nu, _Nu);

				const auto nn = M_k.cols();
				auto Hc_k = Hc.topLeftCorner(nn + _Nu, nn + _Nu);
				Hc_k.topLeftCorner(nn, nn) += M_k.transpose() * Q * M_k;
				Hc_k.topRightCorner(nn, _Nu) += M_k.transpose() * S;
				Hc_k.bottomLeftCorner(_Nu, nn) += ST * M_k;
				Hc_k.bottomRightCorner(_Nu, _Nu) += R;

				auto gc_k = gc.topRows(nn + _Nu);
				gc_k.topRows(nn) += M_k.transpose() * (g_k.topRows(_Nx) + (Q + Q.transpose()) * v);
				gc_k.bottomRows(_Nu) += g_k.bottomRows(_Nu) + (S.transpose() + ST) * v;

				_condensedQP.lb().middleRows(_Nx + k * _Nu, _Nu) = msqp.uMin(k);
				_condensedQP.ub().middleRows(_Nx + k * _Nu, _Nu) = msqp.uMax(k);
			}
			else
			{
				// Final state.
				Hc += M_k.transpose() * H_k * M_k;
				gc += M_k.transpose() * (g_k.topRows(_Nx) + (H_k + H_k.transpose()) * v);
			}
		}
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

		_hotStart = true;

		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw qpOASESException(res);

		/* Get solution of the condensed QP. */
		_problem.getPrimalSolution(_primalCondensedSolution.data());
		//problem.getDualSolution(yOpt);

		// Calculate the solution of the multi-stage QP.
		_primalSolution.topRows(_Nx) = _primalCondensedSolution.topRows(_Nx);
		for (size_type i = 0; i < _Nt; ++i)
		{
			auto z_i = _primalSolution.middleRows(i * _Nz, _Nz);
			auto x_i = z_i.topRows(_Nx);
			auto u_i = z_i.bottomRows(_Nu);
			auto x_i_plus = _primalSolution.middleRows((i + 1) * _Nz, _Nx);
			
			u_i = _primalCondensedSolution.middleRows(_Nx + i * _Nu, _Nu);
			x_i_plus = msqp.C(i) * z_i + msqp.c(i);
		}
	}

	const CondensingSolver::Vector& CondensingSolver::getPrimalSolution() const
	{
		return _primalSolution;
	}

	CondensingSolver::CondensingSolver(size_type nx, size_type nu, size_type nt) 
		: _condensedQP(nx + nt * nu, nt * nx)
		, _Nx(nx), _Nu(nu), _Nt(nt), _Nz(nx + nu)
		, _primalCondensedSolution(nx + nt * nu)
		, _primalSolution(nx + (nx + nu) * nt)
		, _problem(nx + nt * nu, nt * nx)
	{
		qpOASES::Options options;
		_problem.setOptions(options);
	}

	const CondensingSolver::Vector& CondensingSolver::getPrimalCondensedSolution() const
	{
		return _primalCondensedSolution;
	}
}