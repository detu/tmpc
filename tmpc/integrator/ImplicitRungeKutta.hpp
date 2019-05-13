#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/SizeT.hpp>
#include <tmpc/integrator/ButcherTableau.hpp>


namespace tmpc
{
	/**
	 * \defgroup integrators Integrators
	 */

	/**
	 * \brief Implicit Runge-Kutta integrator.
	 * \ingroup integrators
	 *
	 * 
	 */
	template <typename Real>
	class ImplicitRungeKutta
	{
	public:
		ImplicitRungeKutta(size_t nx, size_t nu, ButcherTableau<Real>&& butcher) 
		:	nx_(nx)
		,	nu_(nu)
		,	butcher_(std::move(butcher))
		,	ns_(size(butcher_))
		,	r_(ns_ * nx)
		,	J_(ns_ * nx, ns_ * nx, Real {})
		,	x_(nx)
		,	k_(ns_ * nx_)
		{
		}


		template <typename ODE, typename VT1, typename VT2>
		decltype(auto) operator()(ODE const& ode, Real t0, blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Vector<VT2, blaze::columnVector> const& u, Real h) const
		{
			if (size(x0) != nx_)
				throw std::invalid_argument("Invalid size of x0 in " __FUNCTION__);

			if (size(u) != nu_)
				throw std::invalid_argument("Invalid size of u in " __FUNCTION__);

			// Finding the root of the following equation using Newton method:
			// 0 = f(t_0 + c_i h, x_0 + \sum_{j=1}^s a_{i,j} k_j) - k_i
			k_ = Real {};
			residualMaxNorm_ = inf<Real>();

			for (newtonIterations_ = 0; newtonIterations_ < maxNewtonIterations_; ++newtonIterations_)
			{
				// Calculate the resisuals and the Jacobian
				for (size_t i = 0; i < ns_; ++i)
				{
					x_ = x0;

					for (size_t j = 0; j < ns_; ++j)
						if (a(i, j))
							x_ += a(i, j) * k_(j);

					ode(t0 + c(i) * h, x_, f_, df_dx_);
					r(i) = f_ - k(i);

					for (size_t j = 0; j < ns_; ++j)
						if (a(i, j))
							J(i, j) = a(i, j) * df_dx_;
				}

				residualMaxNorm_ = maxNorm(r_);

				// Solution within tolerance; exit the loop.
				if (residualMaxNorm_ < newtonTolerance_)
					break;

				// Netwon method update: k(n+1) = k(n) - inv(J(n))*r(n)

			}

			noresize(k_[0]) = ode(t0,          x_ = ~x0                , ~u);
			noresize(k_[1]) = ode(t0 + h / 2., x_ = ~x0 + k_[0] * (h / 2.), ~u);
			noresize(k_[2]) = ode(t0 + h / 2., x_ = ~x0 + k_[1] * (h / 2.), ~u);
			noresize(k_[3]) = ode(t0 + h,      x_ = ~x0 + k_[2] * h       , ~u);
	
			return ~x0 + (k_[0] + 2. * k_[1] + 2. * k_[2] + k_[3]) * (h / 6.);
		}

		
		template <typename ODE, typename VT1, bool TF1, typename VT2, bool TF2, 
			typename VT3, bool TF3, typename MT1, bool SO1, typename MT2, bool SO2>
		void operator()(ODE const& ode, Real t0, blaze::Vector<VT1, TF1> const& x0, blaze::Vector<VT2, TF2> const& u, Real h, 
			blaze::Vector<VT3, TF3>& x_next, blaze::Matrix<MT1, SO1>& A, blaze::Matrix<MT2, SO2>& B) const
		{
			if (size(x0) != nx_)
				throw std::invalid_argument("Invalid size of x0 in ImplicitRungeKutta::operator()");

			if (size(u) != nu_)
				throw std::invalid_argument("Invalid size of u in ImplicitRungeKutta::operator()");
				
			throw std::logic_error(__FUNCTION__ " is not implemented");
		}


	private:
		size_t const nx_;
		size_t const nu_;
		ButcherTableau<Real> butcher_;
		size_t const ns_;

		// Vector of residuals
		mutable blaze::DynamicVector<Real, blaze::columnVector> r_;

		// Jacobian
		mutable blaze::DynamicMatrix<Real> J_;

		mutable blaze::DynamicVector<Real, blaze::columnVector> x_;
		mutable blaze::DynamicVector<Real, blaze::columnVector> k_;

		size_t newtonIterations_ = 0;
		size_t maxNewtonIterations_ = 10;
		Real residualMaxNorm_ = inf<Real>();

		// Newton method tolerance
		Real newtonTolerance_ = 1e-10;

		
		Real a(size_t i, size_t j) const
		{
			return butcher_.A()(i, j);
		}


		Real b(size_t j) const
		{
			return butcher_.b()[j];
		}


		Real c(size_t i) const
		{
			return butcher_.c()[i];
		}


		auto k(size_t j)
		{
			return subvector(k_, j * nx_, nx_);
		}


		auto r(size_t j)
		{
			return subvector(r_, j * nx_, nx_);
		}


		auto J(size_t i, size_t j)
		{
			return submatrix(J_, i * nx_, j * nx_, nx_, nx_);
		}
	};
}
