#pragma once

#include <tmpc/integrator/ExplicitIntegrator.hpp>
#include <tmpc/Exception.hpp>


namespace tmpc
{
	/**
	 * \defgroup integrators Integrators
	 */

	/**
	 * \brief Explicit Runge-Kutta integrator.
	 * \ingroup integrators
	 * 
	 */
	template <typename Real>
	class ExplicitRungeKutta
	:	public ExplicitIntegrator<ExplicitRungeKutta<Real>>
	{
	public:
		template <typename Method>
		ExplicitRungeKutta(Method const& method, size_t nx, size_t nu, size_t ny = 0) 
		:	nx_(nx)
		,	nu_(nu)
		,	m_(method.stages())
		,	s_(nx)
		,	S_(nx, nx + nu)
		,	xf_(nx)
		,	Sf_(nx, nx + nu)
		,	k_(m_)
		,	K_(m_)
		,	r_(ny)
		,	J_(ny, nx + nu)
		{
			method.butcherTableau(A_, b_, c_);

			// Check that the Butcher tableau A matrix is strictly lower-triangular
			if (!isStrictlyLower(A_))
				TMPC_THROW_EXCEPTION(std::invalid_argument(
					"Explicit Runge-Kutta methods require Butcher tableau with strictly lower-triangular A matrix"));

			// Resize working variables
			for (size_t i = 0; i < m_; ++i)
			{
				k_[i].resize(nx);
				K_[i].resize(nx, nx + nu);
			}		
		}


		/// @brief Integrate an ODE
		template <typename ODE, typename VT1, typename VT2, typename VT3>
		void operator()(ODE const& ode, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xf) const
		{
			if (size(x0) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			if (size(xf) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of xf"));

			xf_ = ~x0;

			for (size_t i = 0; i < m_; ++i)
			{
				s_ = ~x0;
				for (size_t j = 0; j < i; ++j)
					s_ += h * A_(i, j) * k_[j];

				ode(t0 + h * c_[i], s_, ~u, k_[i]);
				xf_ += h * b_[i] * k_[i];
			}

			~xf = xf_;
		}


		/// @brief Integrate an ODE and calculate forward sensitivities		
		template <typename ODE, 
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2>
		void operator()(ODE const& ode, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Matrix<MT1, SO1> const& S,
			blaze::Vector<VT2, blaze::columnVector> const& u, 
			blaze::Vector<VT3, blaze::columnVector>& xf,
			blaze::Matrix<MT2, SO2>& Sf) const
		{
			if (size(x0) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			if (size(xf) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of xf"));


			xf_ = ~x0;
			Sf_ = ~S;

			for (size_t i = 0; i < m_; ++i)
			{
				s_ = ~x0;
				S_ = ~S;

				for (size_t j = 0; j < i; ++j)
				{
					s_ += h * A_(i, j) * k_[j];
					S_ += h * A_(i, j) * K_[j];
				}

				ode(t0 + h * c_[i], s_, S_, ~u, k_[i], K_[i]);
				xf_ += h * b_[i] * k_[i];
				Sf_ += h * b_[i] * K_[i];
			}

			~xf = xf_;
			~Sf = Sf_;
		}


		/// @brief Integrate an ODE and calculate forward sensitivities,
		/// gradient, and Gauss-Newton Hessian approximation of the integral of least-squares Lagrange term.
		template <typename ODE, 
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2,
			typename VT4,
			typename MT3, bool SO3>
		void operator()(ODE const& ode, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Matrix<MT1, SO1> const& S,
			blaze::Vector<VT2, blaze::columnVector> const& u, 
			blaze::Vector<VT3, blaze::columnVector>& xf,
			blaze::Matrix<MT2, SO2>& Sf,
			Real& l,
			blaze::Vector<VT4, blaze::columnVector>& g,
			blaze::Matrix<MT3, SO3>& H) const
		{
			if (size(x0) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			if (size(xf) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of xf"));


			xf_ = ~x0;
			Sf_ = ~S;

			for (size_t i = 0; i < m_; ++i)
			{
				s_ = ~x0;
				S_ = ~S;

				for (size_t j = 0; j < i; ++j)
				{
					s_ += h * A_(i, j) * k_[j];
					S_ += h * A_(i, j) * K_[j];
				}

				ode(t0 + h * c_[i], s_, S_, ~u, k_[i], K_[i], r_, J_);
				xf_ += h * b_[i] * k_[i];
				Sf_ += h * b_[i] * K_[i];
				
				l += h * b_[i] * sqrNorm(r_) / 2.;
				~g += h * b_[i] * trans(J_) * r_;
				~H += h * b_[i] * trans(J_) * J_;
			}

			~xf = xf_;
			~Sf = Sf_;
		}


	private:
		size_t const nx_;
		size_t const nu_;
		size_t const m_;

		// Butcher Tableau
		blaze::DynamicMatrix<Real> A_;
		blaze::DynamicVector<Real, blaze::rowVector> b_;
		blaze::DynamicVector<Real, blaze::columnVector> c_;
		
		mutable blaze::DynamicVector<Real> s_;
		mutable blaze::DynamicMatrix<Real> S_;

		mutable blaze::DynamicVector<Real> xf_;
		mutable blaze::DynamicMatrix<Real> Sf_;

		mutable std::vector<blaze::DynamicVector<Real>> k_;
		mutable std::vector<blaze::DynamicMatrix<Real>> K_;

		mutable blaze::DynamicVector<Real> r_;
		mutable blaze::DynamicMatrix<Real> J_;
	};
}
