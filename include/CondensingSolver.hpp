#pragma once

#include <MultiStageQP.hpp>
#include <QuadraticProgram.hpp>

namespace camels
{
	class CondensingSolver
	{
	public:
		typedef unsigned size_type;

		CondensingSolver(size_type nx, size_type nu, size_type nt)
			: _qp(nx + nt * nu, nx + nt * nx)
			, _Nx(nx), _Nu(nu), _Nt(nt), _Nz(nx + nu)
		{}

		size_type nIndep() const { return _Nx + _Nu * _Nt; }
		size_type nDep() const { return _Nx * _Nt; }
		size_type nVar() const { return _Nz * _Nt + _Nx; }

		void Condense(const MultiStageQP& msqp)
		{
			assert(msqp.nX() == _Nx && msqp.nU() == _Nu && msqp.nT() == _Nt);
			
			Eigen::MatrixXd M(_Nx, nIndep());
			M.setZero();

			Eigen::VectorXd v(_Nx);
			v.setZero();

			auto& Hc = _qp.H();
			auto& gc = _qp.g();

			Hc.setZero();
			gc.setZero();

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
				}

				_qp.A().middleRows(k * _Nx, _Nx) = M;
				_qp.lbA().middleRows(k * _Nx, _Nx) = msqp.xMin(k) - v;
				_qp.ubA().middleRows(k * _Nx, _Nx) = msqp.xMax(k) - v;

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
				}
				else
				{
					// Final state.
					Hc += M_k.transpose() * H_k * M_k;
					gc += M_k.transpose() * (g_k.topRows(_Nx) + (H_k + H_k.transpose()) * v);
				}
			}
		}

		const QuadraticProgram& QP() const { return _qp; }

	private:
		size_type _Nu;
		size_type _Nx;
		size_type _Nz;
		size_type _Nt;

		QuadraticProgram _qp;
	};
}