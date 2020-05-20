#pragma once

#include <tmpc/integrator/ExplicitIntegrator.hpp>
#include <tmpc/integrator/ImplicitIntegrator.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Test Runge-Kutta integrators using exponential decay as a test system.
    ///
    class DecayTest
	: 	public Test
	{
	protected:
		using Real = double;

		static size_t constexpr NX = 1;
		static size_t constexpr NZ = 0;
		static size_t constexpr NU = 1;
		static size_t constexpr NR = 1;
		

		DecayTest()
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
			
			integrate(~integrator, 
				[this] (auto&&... args) { this->explicitOdeSensitivityResidual(std::forward<decltype(args)>(args)...); },
				0., T_, hMax_, x0_, S_, u_, xf, Sf, l, g, H);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, 1e-10, 1e-3));
			EXPECT_NEAR(l, l_ref_, 1e-6);
			EXPECT_TRUE(approxEqual(g, g_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(H, H_ref_, 1e-10, 1e-4));
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

			(~xdot).resize(NX);
			(~xdot)[0] = -(~u)[0] * (~x)[0];
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

			resize(xdot, NX);
			(~xdot)[0] = -(~u)[0] * (~x)[0];

			resize(Sxdot, NX, NX + NU);
			~Sxdot = {
				{-(~u)[0] * (~Sx)(0, 0), -(~u)[0] * (~Sx)(0, 1) - (~x)[0]}
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

			resize(xdot, NX);
			(~xdot)[0] = -(~u)[0] * (~x)[0];

			resize(Sxdot, NX, NX + NU);
			~Sxdot = {
				{-(~u)[0] * (~Sx)(0, 0), -(~u)[0] * (~Sx)(0, 1) - (~x)[0]}
			};

			~r = ~x;
			~Sr = ~Sx;
		}

	
		/**
		 * \brief Evaluates implicit ODE and Jacobians.
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

			resize(f, NX);
			(~f)[0] = -(~u)[0] * (~x)[0] - (~xdot)[0];

			~Jxdot = -blaze::IdentityMatrix<Real>(NX);

			resize(Jx, NX, NX);
			~Jx = {
				{-(~u)[0]}
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
				{-(~u)[0] * (~Sx)(0, 0), -(~u)[0] * (~Sx)(0, 1) - (~x)[0]}
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

			~r = ~x;
			~Sr = ~Sx;
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

			auto const& u = (~uu)[0];

			// IVP solution
			(~xf).resize(NX);
			(~xf)[0] = exp(-t * u) * (~x0)[0];

			// IVP sensitivities
			(~Sf).resize(NX, NX + NU);
			(~Sf)(0, 0) = exp(-t * u);
			(~Sf)(0, 1) = -exp(-t * u) * t * (~x0)[0];

			// Integral of the Lagrange term
			l = -(((-1. + exp(-2. * t * u)) * pow((~x0)[0], 2)) / (4. * u));

			// Gradient of the integral of the Lagrange term
			resize(g, NX + NU);
			~g = {
				-(((-1. + exp(-2. * t * u)) * (~x0)[0]) / (2. * u)),
				-((exp(-2. * t * u) * (-1. + exp(2. * t * u) - 2. * t * u) * pow((~x0)[0], 2)) / (4. * pow(u, 2)))
			};

			// Integral of the Gauss-Newton Hessian approximation of the Lagrange term
			resize(H, NX + NU, NX + NU);
			(~H)(0, 0) = -((-1. + exp(-2. * t * u)) / (2. * u));
			(~H)(1, 0) = (~H)(0, 1) = -((exp(-2. * t * u) * (-1. + exp(2. * t * u) - 2. * t * u) * (~x0)[0]) / (4. * pow(u, 2)));
			(~H)(1, 1) = ((1. + exp(-2. * t * u) * (-1. - 2. * t * u * (1. + t * u))) * pow((~x0)[0], 2)) / (4. * pow(u, 3));
		}


		void SetUp() override
		{
			analyticalSolution(T_, x0_, u_, xf_ref_, Sf_ref_, l_ref_, g_ref_, H_ref_);
		}


		Real const T_ = 5.;
		Real hMax_ = 0.1;
		blaze::DynamicVector<Real> x0_ {1.};
		blaze::DynamicVector<Real> u_ {1.};

		blaze::DynamicVector<Real> xf_ref_ {NX};
		blaze::DynamicMatrix<Real> Sf_ref_ {NX, NX + NU};
		blaze::DynamicMatrix<Real> S_ {NX, NX + NU, 0.};
		Real l_ref_ {};
		blaze::DynamicVector<Real> g_ref_ {NX + NU};
		blaze::DynamicMatrix<Real> H_ref_ {NX + NU, NX + NU};
	};
}