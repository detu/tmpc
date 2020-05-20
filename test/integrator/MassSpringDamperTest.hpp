#pragma once

#include <tmpc/integrator/ExplicitIntegrator.hpp>
#include <tmpc/integrator/ImplicitIntegrator.hpp>
#include <tmpc/Testing.hpp>

#include <complex>


namespace tmpc :: testing
{
	/// @brief Test Runge-Kutta integrators using mass-spring-damper as a test system.
    ///
    class MassSpringDamperTest
	: 	public Test
	{
	protected:
		using Real = double;
		
		static size_t constexpr NX = 2;
		static size_t constexpr NZ = 0;
		static size_t constexpr NU = 1;
		static size_t constexpr NR = 2;


		MassSpringDamperTest()
		{
			diagonal(S_) = 1.;
		}


		template <typename I>
		void testIntegrate(ExplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf(NX);
			integrate(~integrator, 
				[this] (auto&&... args) { this->explicitOde(std::forward<decltype(args)>(args)...); }, 
				0., T_, hMax_, x0_, u_, xf);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-5));
		}


		template <typename I>
		void testIntegrateWithSensitivities(ExplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf {NX};
			blaze::DynamicMatrix<Real> Sf {NX, NX + NU};
			integrate(~integrator, 
				[this] (auto&&... args) { this->explicitOdeSensitivity(std::forward<decltype(args)>(args)...); }, 
				0., T_, hMax_, x0_, S_, u_, xf, Sf);

			double const tol = 1e-5;
			EXPECT_TRUE(approxEqual(xf, xf_ref_, tol));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, tol));
		}


		template <typename I>
		void testIntegrateLeastSquaresLagrangeTerm(ExplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf(NX);
			blaze::DynamicMatrix<Real> Sf(NX, NX + NU);
			Real l = 0.;
			blaze::DynamicVector<Real> g(NX + NU, 0.);
			blaze::DynamicMatrix<Real> H(NX + NU, NX + NU, 0.);
			
			integrate(integrator, 
				[this] (auto&&... args) { this->explicitOdeSensitivityResidual(std::forward<decltype(args)>(args)...); },
				0., T_, hMax_, x0_, S_, u_, xf, Sf, l, g, H);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, 1e-10, 1e-3));
			EXPECT_NEAR(l, l_ref_, 1e-6);
			EXPECT_TRUE(approxEqual(g, g_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(H, H_ref_, 1e-10, 1e-5));
		}


		template <typename I>
		void testIntegrate(ImplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf(NX);
			integrate(~integrator,
				[this] (auto&&... args) { this->implicitOde(std::forward<decltype(args)>(args)...); },
				0., T_, hMax_, x0_, u_, xf);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-5));
		}


		template <typename I>
		void testIntegrateWithSensitivities(ImplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf {NX};
			blaze::DynamicMatrix<Real> Sf {NX, NX + NU};
			integrate(~integrator, 
				[this] (auto&&... args) { this->implicitOde(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->implicitOdeSensitivity(std::forward<decltype(args)>(args)...); },
				0., T_, hMax_, x0_, S_, u_, xf, Sf);

			double const tol = 1e-5;
			EXPECT_TRUE(approxEqual(xf, xf_ref_, tol));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, tol));
		}


		template <typename I>
		void testIntegrateLeastSquaresLagrangeTerm(ImplicitIntegrator<I> const& integrator)
		{
			blaze::DynamicVector<Real> xf(NX);
			blaze::DynamicMatrix<Real> Sf(NX, NX + NU);
			Real l = 0.;
			blaze::DynamicVector<Real> g(NX + NU, 0.);
			blaze::DynamicMatrix<Real> H(NX + NU, NX + NU, 0.);
			
			integrate(~integrator, 
				[this] (auto&&... args) { this->implicitOde(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->implicitOdeSensitivity(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->residual(std::forward<decltype(args)>(args)...); },
				0., T_, hMax_, x0_, S_, u_, xf, Sf, l, g, H);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, 1e-10, 1e-3));
			EXPECT_NEAR(l, l_ref_, 1e-6);
			EXPECT_TRUE(approxEqual(g, g_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(H, H_ref_, 1e-10, 1e-5));
		}
		

	private:
		/**
		 * \brief Evaluates ODE.
		 */
		template <typename VT1,	typename VT2, typename VT3>
		void explicitOde(double t, 
			blaze::Vector<VT1, blaze::columnVector> const& x,
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xdot) const
		{
			if (size(x) != NX || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			auto const& k = k_;
			auto const& kappa = kappa_;

			(~xdot).resize(NX);
			(~xdot)[0] = (~x)[1];
			(~xdot)[1] = -k * (~x)[0] - 2. * kappa * (~x)[1] + (~u)[0];
		}


		/**
		 * \brief Evaluates ODE and sensitivities.
		 */
		template <
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2>
		void explicitOdeSensitivity(double t, 
			blaze::Vector<VT1, blaze::columnVector> const& x,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xdot,
			blaze::Matrix<MT2, SO2>& Sxdot) const
		{
			if (size(x) != NX || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			auto const& k = k_;
			auto const& kappa = kappa_;

			~xdot = {
				(~x)[1],
				-k * (~x)[0] - 2. * kappa * (~x)[1] + (~u)[0]
			};

			~Sxdot = {
				{(~Sx)(1, 0), (~Sx)(1, 1), (~Sx)(1, 2)},
				{-k * (~Sx)(0, 0) - 2. * kappa * (~Sx)(1, 0),
				-k * (~Sx)(0, 1) - 2. * kappa * (~Sx)(1, 1),
				1. - k * (~Sx)(0, 2) - 2. * kappa * (~Sx)(1, 2)}
			};
		}


		/**
		 * \brief Evaluates ODE, ODE sensitivities, resudual, and Jacobian of the residual.
		 */
		template <
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename VT3,
			typename MT2, bool SO2,
			typename VT4,
			typename MT3, bool SO3>
		void explicitOdeSensitivityResidual(double t, 
			blaze::Vector<VT1, blaze::columnVector> const& x,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xdot,
			blaze::Matrix<MT2, SO2>& Sxdot,
			blaze::Vector<VT4, blaze::columnVector>& r,
			blaze::Matrix<MT3, SO3>& Sr) const
		{
			if (size(x) != NX || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			auto const& k = k_;
			auto const& kappa = kappa_;

			~xdot  = {
				(~x)[1],
				-k * (~x)[0] - 2. * kappa * (~x)[1] + (~u)[0]
			};

			~Sxdot = {
				{(~Sx)(1, 0), (~Sx)(1, 1), (~Sx)(1, 2)},
				{-k * (~Sx)(0, 0) - 2. * kappa * (~Sx)(1, 0),
				-k * (~Sx)(0, 1) - 2. * kappa * (~Sx)(1, 1),
				1. - k * (~Sx)(0, 2) - 2. * kappa * (~Sx)(1, 2)}
			};

			~r = {
				(~x)[1] * (~u)[0],
				(~u)[0]
			};

			~Sr = {
				{(~u)[0] * (~Sx)(1, 0), (~u)[0] * (~Sx)(1, 1), (~u)[0] * (~Sx)(1, 2) + (~x)[1]},
				{0., 0., 1.}
			};
		}


		/**
		 * \brief Evaluates the implicit ODE and its Jacobians.
		 */
		template <
			typename VT1,
			typename VT2,
			typename VT3,
			typename VT4,
			typename VT5,
			typename MT1, bool SO1,
			typename MT2, bool SO2,
			typename MT3, bool SO3>
		void implicitOde(
			Real t,
			blaze::Vector<VT1, blaze::columnVector> const& xdot,
			blaze::Vector<VT2, blaze::columnVector> const& x,
			blaze::Vector<VT3, blaze::columnVector> const& z,
			blaze::Vector<VT4, blaze::columnVector> const& u,
			blaze::Vector<VT5, blaze::columnVector>& f,
			blaze::Matrix<MT1, SO1>& Jxdot,
			blaze::Matrix<MT2, SO2>& Jx,
			blaze::Matrix<MT3, SO3>& Jz) const
		{
			if (size(xdot) != NX || size(x) != NX || size(z) != NZ || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			auto const& k = k_;
			auto const& kappa = kappa_;

			resize(f, NX);
			(~f)[0] = (~x)[1] - (~xdot)[0];
			(~f)[1] = -k * (~x)[0] - 2. * kappa * (~x)[1] + (~u)[0] - (~xdot)[1];

			~Jxdot = -blaze::IdentityMatrix<Real>(NX);

			resize(Jx, NX, NX);
			~Jx = {
				{0., 1.},
				{-k, -2. * kappa}
			};
		}


		/**
		 * \brief Evaluates the implicit ODE sensitivity.
		 */
		template <
			typename VT1,
			typename VT2,
			typename MT1, bool SO1,
			typename VT3,
			typename VT4,
			typename MT2, bool SO2>
		void implicitOdeSensitivity(double t,
			blaze::Vector<VT1, blaze::columnVector> const& xdot,
			blaze::Vector<VT2, blaze::columnVector> const& x,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT3, blaze::columnVector> const& z,
			blaze::Vector<VT4, blaze::columnVector> const& u,
			blaze::Matrix<MT2, SO2>& Sf) const
		{
			if (size(xdot) != NX || size(x) != NX || size(z) != NZ || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			resize(Sf, NX, NX + NU);
			~Sf = {
				{(~Sx)(1, 0), (~Sx)(1, 1), (~Sx)(1, 2)},
				{-k_ * (~Sx)(0, 0) - 2. * kappa_ * (~Sx)(1, 0),
				-k_ * (~Sx)(0, 1) - 2. * kappa_ * (~Sx)(1, 1),
				1. - k_ * (~Sx)(0, 2) - 2. * kappa_ * (~Sx)(1, 2)}
			};
		}
			

		/**
		 * \brief Evaluates the residual and Jacobian of the residual.
		 */
		template <
			typename VT1,
			typename MT1, bool SO1,
			typename VT2,
			typename MT2, bool SO2,
			typename VT3,
			typename VT4,
			typename MT3, bool SO3
		>
		void residual(double t, 
			blaze::Vector<VT1, blaze::columnVector> const& x,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT2, blaze::columnVector> const& z,
			blaze::Matrix<MT2, SO2> const& Sz,
			blaze::Vector<VT3, blaze::columnVector> const& u,
			blaze::Vector<VT4, blaze::columnVector>& r,
			blaze::Matrix<MT3, SO3>& Sr) const
		{
			if (size(x) != NX || size(z) != NZ || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			~r = {
				(~x)[1] * (~u)[0],
				(~u)[0]
			};

			~Sr = {
				{(~u)[0] * (~Sx)(1, 0), (~u)[0] * (~Sx)(1, 1), (~u)[0] * (~Sx)(1, 2) + (~x)[1]},
				{0., 0., 1.}
			};
		}


		/**
		 * \brief Calculate analytical solution of the ODE and its analytical sensitivities.
		 */
		template <
			typename VT1,
			typename VT2,
			typename VT3,
			typename MT1, bool SO1, 
			typename VT4,
			typename MT2, bool SO2>
		void analyticalSolution(double t, 
			blaze::Vector<VT1, blaze::columnVector> const& x0,
			blaze::Vector<VT2, blaze::columnVector> const& uu,
			blaze::Vector<VT3, blaze::columnVector>& xf,
			blaze::Matrix<MT1, SO1>& Sf,
			double& l,
			blaze::Vector<VT4, blaze::columnVector>& g,
			blaze::Matrix<MT2, SO2>& H) const
		{
			using std::pow;

			if (size(x0) != NX || size(uu) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			auto const& q0 = (~x0)[0];
			auto const& v0 = (~x0)[1];
			auto const& u = (~uu)[0];

			auto const s = sqrt(std::complex<double> {-k_ + pow(kappa_, 2)});
			auto const qf = (1./(2.*k_*s)) * exp(-t*(s+kappa_)) 
				* (-u*((1.+exp(2.*s*t)-2.*exp(t*(s+kappa_)))*s+(-1.+exp(2.*s*t))*kappa_)
				+k_*((-1.+exp(2.*s*t))*v0-q0*kappa_+q0*(s+exp(2.*s*t)*(s+kappa_))));
			auto const vf = (exp(-t*kappa_)*(s*v0*cosh(s*t) + (-k_*q0 + u - v0*kappa_)*sinh(s*t)))/s;

			// IVP solution
			(~xf).resize(NX);
			(~xf)[0] = real(qf);
			(~xf)[1] = real(vf);

			// IVP sensitivities
			(~Sf).resize(NX, NX + NU);
			(~Sf)(0, 0) = real((exp(-t*kappa_)*(s*cosh(s*t) + kappa_*sinh(s*t)))/s);
			(~Sf)(0, 1) = real((exp(-t*(s + kappa_))*(-1. + exp(2.*s*t)))/(2.*s));
			(~Sf)(0, 2) = real((s - exp(-t*kappa_)*(s*cosh(s*t) + kappa_*sinh(s*t)))/(k_*s));
			(~Sf)(1, 0) = real(-((exp(-t*kappa_)*k_*sinh(s*t))/s));
			(~Sf)(1, 1) = real((exp(-t*kappa_)*(s*cosh(s*t) - kappa_*sinh(s*t)))/s);
			(~Sf)(1, 2) = real((exp(-t*kappa_)*sinh(s*t))/s);

			// Lagrange term integral
			l = real((1./(8.*pow(s, 2)*(s-kappa_)*kappa_*(s+kappa_)))*exp(-2.*t*kappa_)*pow(u, 2)
				*((s-kappa_)*(s+kappa_)*(k_*q0-u+v0*(-s+kappa_))*(k_*q0-u+v0*(s+kappa_))-exp(2.*t*kappa_)*pow(s, 2)
				*(pow(-k_*q0+u, 2)+(pow(v0, 2)+4.*t*kappa_)*(-pow(s, 2)+pow(kappa_, 2)))+kappa_*(pow(k_, 2)*pow(q0, 2)*kappa_
				+kappa_*pow(u-v0*kappa_, 2)+pow(s, 2)*v0*(2.*u-v0*kappa_)-2.*k_*q0*(pow(s, 2)*v0+kappa_*(u-v0*kappa_)))*cosh(2.*s*t)
				+s*kappa_*(pow(-k_*q0+u, 2)+pow(v0, 2)*(s-kappa_)*(s+kappa_))*sinh(2.*s*t)));

			// Lagrange term integral gradient
			~g = {
				real((1./(4.*pow(s, 2)*kappa_*(-s + kappa_)*(s + kappa_)))*exp(-2.*t*kappa_)
					*k_*pow(u, 2)*(exp(2.*t*kappa_)*pow(s, 2)*(k_*q0 - u) + (k_*q0 - u + v0*kappa_)*(-pow(s, 2) + pow(kappa_, 2)) 
					+ kappa_*(pow(s, 2)*v0 + kappa_*(-k_*q0 + u - v0*kappa_))*cosh(2.*s*t) + s*(-k_*q0 + u)*kappa_*sinh(2.*s*t))),
				real((1./(4.*pow(s, 2)*kappa_))*exp(-2.*t*kappa_)*pow(u, 2)*((-1. + exp(2.*t*kappa_))*pow(s, 2)*v0 + kappa_*(k_*q0 - u
					+ v0*kappa_) + kappa_*(-k_*q0 + u - v0*kappa_)*cosh(2.*s*t) + s*v0*kappa_*sinh(2.*s*t))),
				real((1./(4.*pow(s, 2)*kappa_*(-s + kappa_)*(s + kappa_)))*exp(-2.*t*kappa_)*u*((s - kappa_)*(s + kappa_)*(pow(s, 2)*pow(v0, 2) 
					- (k_*q0 - 2.*u + v0*kappa_)*(k_*q0 - u + v0*kappa_)) + exp(2.*t*kappa_)*pow(s, 2)*((k_*q0 - 2.*u)*(k_*q0 - u) + (pow(v0, 2) + 
					4.*t*kappa_)*(-pow(s, 2) + pow(kappa_, 2))) - kappa_*(pow(k_, 2)*pow(q0, 2)*kappa_ + pow(s, 2)*v0*(3.*u - v0*kappa_)
					+ kappa_*(-2.*u + v0*kappa_)*(-u + v0*kappa_) + k_*q0*(-2.*pow(s, 2)*v0 + kappa_*(-3.*u + 2.*v0*kappa_)))*cosh(2.*s*t) - 
					s*kappa_*((k_*q0 - 2.*u)*(k_*q0 - u) + pow(v0, 2)*(s - kappa_)*(s + kappa_))*sinh(2.*s*t)))
			};

			// Lagrange term Gauss-Newton Hessian approximation
			(~H)(0, 0) = real((exp(-2.*t*kappa_)*pow(k_, 2)*pow(u, 2)*((-1. + exp(
				2.*t*kappa_))*pow(s, 2) + pow(kappa_, 2) - kappa_*(kappa_*cosh(2.*s*t) + s*sinh(2.*s*t))))
				/(4.*pow(s, 2)*(-pow(s, 2)*kappa_ + pow(kappa_, 3))));
			(~H)(1, 0) = (~H)(0, 1) = real(-((exp(-2.*t*kappa_)*k_*pow(u, 2)*pow(sinh(s*t), 2))/(2.*pow(s, 2))));
			(~H)(1, 1) = real((exp(-2.*t*kappa_)*pow(u, 2)*((-1. + exp(
				2.*t*kappa_))*pow(s, 2) + pow(kappa_, 2) - pow(kappa_, 2)*cosh(2.*s*t) + 
				s*kappa_*sinh(2.*s*t)))/(4.*pow(s, 2)*kappa_));
			(~H)(2, 0) = (~H)(0, 2) = real((exp(-2.*t*kappa_)*k_*u*(exp(2.*t*kappa_)
				*pow(s, 2)*(k_*q0 - 2.*u) + (k_*q0 - 2.*u + 
				v0*kappa_)*(-pow(s, 2) + pow(kappa_, 2)) + kappa_*(pow(s, 2)*v0 - kappa_
				*(k_*q0 - 2.*u + v0*kappa_))*cosh(2.*s*t) + 
				s*(-k_*q0 + 2.*u)*kappa_*sinh(2.*s*t)))/(4.*pow(s, 2)*kappa_*(-s + kappa_)*(s + kappa_)));
			(~H)(2, 1) = (~H)(1, 2) = real((exp(-2.*t*kappa_)*u*((-1. + exp(
				2.*t*kappa_))*pow(s, 2)*v0 + kappa_*(k_*q0 - 2.*u + 
				v0*kappa_) - kappa_*(k_*q0 - 2.*u + v0*kappa_)*cosh(2.*s*t) +
				s*v0*kappa_*sinh(2.*s*t)))/(4.*pow(s, 2)*kappa_));
			(~H)(2, 2) = real((exp(-2.*t*kappa_)*((s - kappa_)*(s + kappa_)*(k_*q0 - 2.*u + 
				v0*(-s + kappa_))*(k_*q0 - 2.*u + v0*(s + kappa_)) + 
				exp(2.*t*kappa_)*pow(s, 2)*(-pow(k_*q0 - 2.*u, 2) + (s - kappa_)*(s + kappa_)*(pow(v0, 2) + 
				4.*t*kappa_)) + kappa_*(pow(k_, 2)*pow(q0, 2)*kappa_ + 
				pow(s, 2)*v0*(4.*u - v0*kappa_) + kappa_*pow(-2.*u + v0*kappa_, 2) - 
				2.*k_*q0*(2.*u*kappa_ + v0*(s - kappa_)*(s + kappa_)))*cosh(2.*s*t) + 
				s*kappa_*(pow(k_*q0 - 2.*u, 2) + 
				pow(v0, 2)*(s - kappa_)*(s + kappa_))*sinh(2.*s*t)))/(4.*pow(s, 2)*(s - kappa_)*kappa_*(s + kappa_)));
		}


		void SetUp() override
		{
			analyticalSolution(T_, x0_, u_, xf_ref_, Sf_ref_, l_ref_, g_ref_, H_ref_);
		}


		Real k_ = 1.;
		Real kappa_ = 0.2;
		Real const T_ = 5.;
		Real hMax_ = 0.1;
		blaze::DynamicVector<Real> x0_ {1., 1.};
		blaze::DynamicVector<Real> u_ {1.};

		blaze::DynamicVector<Real> xf_ref_ {NX};
		blaze::DynamicMatrix<Real> Sf_ref_ {NX, NX + NU};
		blaze::DynamicMatrix<Real> S_ {NX, NX + NU, 0.};
		Real l_ref_ {};
		blaze::DynamicVector<Real> g_ref_ {NX + NU};
		blaze::DynamicMatrix<Real> H_ref_ {NX + NU, NX + NU};
	};
}