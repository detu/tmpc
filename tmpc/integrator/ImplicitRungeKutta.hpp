#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/SizeT.hpp>
#include <tmpc/integrator/ButcherTableau.hpp>
#include <tmpc/numeric/NewtonSolver.hpp>

#include <boost/throw_exception.hpp>

#include <stdexcept>


namespace tmpc
{
	/**
	 * \defgroup integrators Integrators
	 */

	/**
	 * \brief Implicit Runge-Kutta integrator.
	 * \ingroup integrators
	 * 
	 */
	template <typename Real>
	class ImplicitRungeKutta
	{
	public:
		ImplicitRungeKutta(size_t nx, size_t nu, ButcherTableau<Real> butcher) 
		:	nx_(nx)
		,	nu_(nu)
		,	butcher_(std::move(butcher))
		,	ns_(size(butcher_))
		,	r_(ns_ * nx)
		,	x_(nx)
		,	f_(nx)
		,	df_dx_(nx, nx)
		,	k_(ns_ * nx_)
		,	newtonSolver_(ns_ * nx)
		{
		}


		template <typename ODE, typename VT1, typename VT2>
		auto const& operator()(ODE const& ode, Real t0, blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Vector<VT2, blaze::columnVector> const& u, Real h) const
		{
			return (*this)(ode, t0, x0, u, h, [] (size_t, auto const&, auto const&, auto const&) {});
		}


		template <typename ODE, typename VT1, typename VT2, typename Monitor>
		auto const& operator()(ODE const& ode, Real t0, blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Vector<VT2, blaze::columnVector> const& u, Real h, Monitor monitor) const
		{
			if (size(x0) != nx_)
				BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			// Finding the root of the following equation using Newton method:
			// 0 = f(t_0 + c_i h, x_0 + \sum_{j=1}^s a_{i,j} k_j) - k_i
			k_ = newtonSolver_.solve(
				[&] (auto const& k, auto& r, auto& J)
				{
					for (size_t i = 0; i < ns_; ++i)
					{
						x_ = Real {};
						for (size_t j = 0; j < ns_; ++j)
							x_ += a(i, j) * subvector(k, j * nx_, nx_);
						x_ = h * x_ + ~x0;

						ode(t0 + c(i) * h, x_, ~u, f_, df_dx_);
						subvector(r, i * nx_, nx_) = f_ - subvector(k, i * nx_, nx_);

						for (size_t j = 0; j < ns_; ++j)
							submatrix(J, i * nx_, j * nx_, nx_, nx_) = (h * a(i, j)) * df_dx_;
					}

					J -= blaze::IdentityMatrix<Real>(ns_ * nx_);
				}, 
				blaze::ZeroVector<Real>(ns_ * nx_),
				monitor
			);

			// Calculating the value of the integral
			x_ = Real {};
			for (size_t i = 0; i < ns_; ++i)
				x_ += b(i) * subvector(k_, i * nx_, nx_);

			x_ *= h;
			x_ += x0;
	
			return x_;
		}

		
		template <typename ODE, typename VT1, bool TF1, typename VT2, bool TF2, 
			typename VT3, bool TF3, typename MT1, bool SO1, typename MT2, bool SO2>
		void operator()(ODE const& ode, Real t0, blaze::Vector<VT1, TF1> const& x0, blaze::Vector<VT2, TF2> const& u, Real h, 
			blaze::Vector<VT3, TF3>& x_next, blaze::Matrix<MT1, SO1>& A, blaze::Matrix<MT2, SO2>& B) const
		{
			if (size(x0) != nx_)
				BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));
				
			throw BOOST_THROW_EXCEPTION(std::logic_error("Function is not implemented"));
		}


		/// @brief Get max number of Newton iterations
		size_t newtonMaxIterations() const
		{
			return newtonSolver_.maxIterations();
		}


		/// @brief Set max number of Newton iterations
		void newtonMaxIterations(size_t val)
		{
			newtonSolver_.maxIterations(val);
		}


		/// @brief Set Newton method backtracking parameter
		void newtonBacktrackingAlpha(Real val)
		{
			newtonSolver_.backtrackingAlpha(val);
		}


		/// @brief Get number of Newton iterations made on the last operator() call.
		size_t newtonIterations() const
		{
			return newtonSolver_.iterations();
		}


	private:
		size_t const nx_;
		size_t const nu_;
		ButcherTableau<Real> butcher_;
		size_t const ns_;

		// Vector of residuals
		mutable blaze::DynamicVector<Real, blaze::columnVector> r_;

		mutable blaze::DynamicVector<Real, blaze::columnVector> x_;
		mutable blaze::DynamicVector<Real, blaze::columnVector> f_;
		mutable blaze::DynamicMatrix<Real> df_dx_;
		mutable blaze::DynamicVector<Real, blaze::columnVector> k_;

		mutable NewtonSolver<Real> newtonSolver_;


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
	};
}
