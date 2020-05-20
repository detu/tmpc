#pragma once

#include <tmpc/integrator/ImplicitIntegrator.hpp>
#include <tmpc/numeric/NewtonSolver.hpp>

#include <tmpc/Exception.hpp>

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
	:	public ImplicitIntegrator<ImplicitRungeKutta<Real>>
	{
	public:
		template <typename Method>
		ImplicitRungeKutta(Method const& method, size_t nx, size_t nz, size_t nu, size_t ny = 0)
		:	nx_(nx)
		,	nz_(nz)
		,	nw_(nx + nz)
		,	nu_(nu)
		,	ny_(ny)
		,	m_(method.stages())
		,	s_(m_ * nx)
		,	S_(m_ * nx, nx + nu)
		,	df_dx_(nw_, nx)
		,	kz_(m_ * nw_, Real {})
		,	K_(m_ * nw_, nx + nu)
		,	r_(ny)
		,	Jr_(ny, nx + nu)
		,	newtonSolver_(m_ * nw_)
		{
			method.butcherTableau(A_, b_, c_);
		}


		template <typename DAE, typename VT1, typename VT2, typename VT3>
		void operator()(DAE const& dae, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xf) const
		{
			(*this)(dae, t0, h, ~x0, ~u, ~xf, [] (size_t, auto const&, auto const&, auto const&) {});
		}


		template <typename DAE, typename VT1, typename VT2, typename VT3, typename Monitor>
		void operator()(DAE const& dae, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0, 
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xf,
			Monitor monitor) const
		{
			if (size(x0) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			// Finding the root of the following equation using Newton method:
			// 0 = f(t_0 + c_i h, x_0 + \sum_{j=1}^s a_{i,j} k_j) - k_i

			// If warm-starting, reuse the previous value of k as the starting point,
			// otherwise reset it to 0.
			if (!warmStart_)
				reset(kz_);

			newtonSolver_(newtonResidual(dae, t0, h, ~x0, ~u), kz_, kz_, monitor);

			// Calculating the value of the integral
			~xf = x0;
			for (size_t i = 0; i < m_; ++i)
				~xf += h * b_[i] * subvector(kz_, i * nw_, nx_);
		}

		
		template <
			typename DAE, 
			typename DAE_S,
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2>
		void operator()(
			DAE const& dae, 
			DAE_S const& dae_s,
			Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xf,
			blaze::Matrix<MT2, SO2>& Sf) const
		{
			if (size(x0) != nx_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of x0"));

			if (size(u) != nu_)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of u"));

			// Finding the root of the following equation using Newton method:
			// 0 = f(t_0 + c_i h, x_0 + \sum_{j=1}^s a_{i,j} k_j) - k_i

			// If warm-starting, reuse the previous value of k as the starting point,
			// otherwise reset it to 0.
			if (!warmStart_)
				reset(kz_);

			// Calculate implicit DAE/DAE solution k and its sensitivities
			newtonSolver_(
				newtonResidual(dae, t0, h, ~x0, ~u),
				newtonParamSensitivity(dae_s, t0, h, ~x0, ~Sx, ~u),
				kz_, kz_, K_);

			// Calculating sensitivities of intermediate state variables.
			for (size_t i = 0; i < m_; ++i)
			{
				S(i) = ~Sx;
				
				for (size_t j = 0; j < m_; ++j)
					S(i) += h * A_(i, j) * submatrix(K_, j * nw_, 0, nx_, nx_ + nu_);
			}

			// Calculating the value of the integral and final state sensitivities
			~xf = x0;
			~Sf = ~Sx;

			for (size_t i = 0; i < m_; ++i)
			{
				~xf += h * b_[i] * subvector(kz_, i * nw_, nx_);
				~Sf += h * b_[i] * submatrix(K_, i * nw_, 0, nx_, nx_ + nu_);
			}
		}


		template <
			typename DAE,
			typename DAE_S,
			typename Residual,
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2,
			typename VT4,
			typename MT3, bool SO3>
		void operator()(
			DAE const& dae,
			DAE_S const& dae_s,
			Residual const& res,
			Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Matrix<MT1, SO1> const& Sx,
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

			// Calculate the final state value and its sensitivities
			(*this)(dae, dae_s, t0, h, ~x0, ~Sx, ~u, ~xf, ~Sf);

			// Calculate the integral of the Lagrange term, its gradient, and its Gauss-Newton Hessian
			for (size_t i = 0; i < m_; ++i)
			{
				// Calculate the residual and its sensitivities at i-th intermediate point
				res(t0 + h * c_[i], s(i), S(i), z(i), Z(i), ~u, r_, Jr_);

				// Update the cost, its gradient, and its Gauss-Newton Hessian
				l += h * b_[i] * sqrNorm(r_) / 2.;
				~g += h * b_[i] * trans(Jr_) * r_;
				~H += h * b_[i] * trans(Jr_) * Jr_;
			}
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


		/// @brief Get warm start value.
		bool warmStart() const
		{
			return warmStart_;
		}


		/// @brief Switch warm start on and off.
		void warmStart(bool val)
		{
			warmStart_ = val;
		}


	private:
		size_t const nx_;
		size_t const nz_;
		size_t const nw_;
		size_t const nu_;
		size_t const ny_;
		size_t const m_;

		// Butcher Tableau
		blaze::DynamicMatrix<Real> A_;
		blaze::DynamicVector<Real, blaze::rowVector> b_;
		blaze::DynamicVector<Real, blaze::columnVector> c_;

		// Intermediate state variables
		mutable blaze::DynamicVector<Real, blaze::columnVector> s_;
		
		// Sensitivity of s w.r.t. (x,u)
		mutable blaze::DynamicMatrix<Real> S_;
		
		mutable blaze::DynamicMatrix<Real, blaze::columnMajor> df_dx_;
		mutable blaze::DynamicVector<Real, blaze::columnVector> kz_;

		// Sensitivity of K w.r.t. (x,u)
		mutable blaze::DynamicMatrix<Real, blaze::columnMajor> K_;

		// Holds the residual
		mutable blaze::DynamicVector<Real> r_;

		// Holds the Jacobian of the residual
		mutable blaze::DynamicMatrix<Real> Jr_;

		mutable NewtonSolver<Real> newtonSolver_;

		// Use previous solution as the initial point in the Newton method
		bool warmStart_ = false;


		template <typename DAE, typename VT1, typename VT2>
		auto newtonResidual(DAE const& dae, Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Vector<VT2, blaze::columnVector> const& u) const
		{
			return [this, &dae, t0, h, &x0, &u] (auto const& kz, auto& r, auto& J)
			{
				for (size_t i = 0; i < m_; ++i)
				{
					s(i) = ~x0;
					for (size_t j = 0; j < m_; ++j)
						s(i) += h * A_(i, j) * subvector(kz, j * nw_, nx_);

					auto const k_i = subvector(kz, i * nw_, nx_);
					auto const z_i = subvector(kz, i * nw_ + nx_, nz_);
					auto f = subvector(r, i * nw_, nw_);
					auto Jxdot = submatrix(J, i * nw_, i * nw_, nw_, nx_);
					auto Jz = submatrix(J, i * nw_, i * nw_ + nx_, nw_, nz_);
					dae(t0 + c_[i] * h, k_i, s(i), z_i, ~u, f, Jxdot, df_dx_, Jz);

					for (size_t j = 0; j < m_; ++j)
					{
						auto Jx_ij = submatrix(J, i * nw_, j * nw_, nw_, nx_);
						auto Jz_ij = submatrix(J, i * nw_, j * nw_ + nx_, nw_, nz_);

						if (j != i)
						{
							Jx_ij = h * A_(i, j) * df_dx_;
							reset(Jz_ij);
						}
						else
						{
							Jx_ij += h * A_(i, j) * df_dx_;
						}						
					}
				}
			};
		}


		template <
			typename DAE_S,
			typename VT1,
			typename MT, bool TF,
			typename VT2
		>
		auto newtonParamSensitivity(
			DAE_S const& dae_s,
			Real t0, Real h,
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Matrix<MT, TF> const& Sx,
			blaze::Vector<VT2, blaze::columnVector> const& u) const
		{
			return [this, &dae_s, t0, h, &Sx, &u] (auto const& kz, auto& df_dp)
			{
				for (size_t i = 0; i < m_; ++i)
				{
					// NOTE: s(i) have valid values here,
					// which have already been calculated while calculating the Newton residual.
					auto const k_i = subvector(kz, i * nw_, nx_);
					auto const z_i = subvector(kz, i * nw_ + nx_, nz_);
					auto dfi_dxu = submatrix(df_dp, i * nw_, 0, nw_, nx_ + nu_);
					dae_s(t0 + c_[i] * h, k_i, s(i), ~Sx, z_i, ~u, dfi_dxu);
				}
			};
		}


		decltype(auto) s(size_t i) const
		{
			return subvector(s_, nx_ * i, nx_);
		}


		decltype(auto) S(size_t i) const
		{
			return submatrix(S_, nx_ * i, 0, nx_, nx_ + nu_);
		}


		decltype(auto) z(size_t i) const
		{
			return subvector(kz_, nw_ * i + nx_, nz_);
		}


		decltype(auto) Z(size_t i) const
		{
			return submatrix(K_, nw_ * i + nx_, 0, nz_, nx_ + nu_);
		}
	};
}
