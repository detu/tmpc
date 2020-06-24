#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Test integration of a simple DAE of the following form:
	///
	/// 0 = \dot{x} + z
	/// 0 = x - z^2
	///
	class IrkSimpleDaeTest
	:	public Test
	{
	protected:
		using Real = double;
		
		static size_t constexpr NX = 1;
		static size_t constexpr NZ = 1;
		static size_t constexpr NU = 0;
		static size_t constexpr NR = 2;

		
		IrkSimpleDaeTest()
		{
			diagonal(S_) = 1.;
		}
		

		template <typename Method>
		void testIntegrate(Method const& method, double abs_tol, double rel_tol)
		{
			ImplicitRungeKutta<Real> irk(method, NX, NZ, NU);
			irk.warmStart(true);

			VecX xf;
			integrate(irk,
				[this] (auto&&... args) { this->implicitDae(std::forward<decltype(args)>(args)...); },
				0., T_, numSteps_, x0_, u_, xf
			);

			TMPC_EXPECT_APPROX_EQ(xf, xf_ref_, abs_tol, rel_tol);
		}


		template <typename Method>
		void testIntegrateWithSensitivities(Method const& method, double abs_tol, double rel_tol)
		{
			ImplicitRungeKutta<Real> irk(method, NX, NZ, NU);
			irk.warmStart(true);

			VecX xf;
			blaze::DynamicMatrix<Real> Sf(NX, NX + NU);
			integrate(irk,
				[this] (auto&&... args) { this->implicitDae(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->implicitDaeSensitivity(std::forward<decltype(args)>(args)...); },
				0., T_, numSteps_, x0_, S_, u_, xf, Sf
			);

			TMPC_EXPECT_APPROX_EQ(xf, xf_ref_, abs_tol, rel_tol);
			TMPC_EXPECT_APPROX_EQ(Sf, Sf_ref_, abs_tol, rel_tol);
		}


		template <typename Method>
		void testIntegrateLeastSquaresLagrangeTerm(Method const& method, double abs_tol, double rel_tol)
		{
			ImplicitRungeKutta<Real> irk(method, NX, NZ, NU, NR);
			irk.warmStart(true);

			blaze::DynamicVector<Real> xf(NX);
			blaze::DynamicMatrix<Real> Sf(NX, NX + NU);
			Real l = 0.;
			blaze::DynamicVector<Real> g(NX + NU, 0.);
			blaze::DynamicMatrix<Real> H(NX + NU, NX + NU, 0.);
			
			integrate(irk, 
				[this] (auto&&... args) { this->implicitDae(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->implicitDaeSensitivity(std::forward<decltype(args)>(args)...); },
				[this] (auto&&... args) { this->residual(std::forward<decltype(args)>(args)...); },
				0., T_, numSteps_, x0_, S_, u_, xf, Sf, l, g, H);

			EXPECT_TRUE(approxEqual(xf, xf_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(Sf, Sf_ref_, 1e-10, 1e-3));
			EXPECT_NEAR(l, l_ref_, 1e-6);
			EXPECT_TRUE(approxEqual(g, g_ref_, 1e-10, 1e-5));
			EXPECT_TRUE(approxEqual(H, H_ref_, 1e-10, 1e-5));
		}


	private:
		using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
		using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;

		VecX const x0_ {1.};
		VecU const u_ {};
		Real const T_ = 1.;
		size_t const numSteps_ = 100;
		
		VecX xf_ref_;
		blaze::DynamicMatrix<Real> Sf_ref_ {NX, NX + NU};
		blaze::DynamicMatrix<Real> S_ {NX, NX + NU, 0.};
		Real l_ref_ {};
		blaze::DynamicVector<Real> g_ref_ {NX + NU};
		blaze::DynamicMatrix<Real> H_ref_ {NX + NU, NX + NU};


		void SetUp() override
		{
			// Calculate the correct values analytically
			xf_ref_ = pow(2. * sqrt(x0_) - T_, 2) / 4.;
			Sf_ref_ = {{1. - T_ / (2. * sqrt(x0_[0]))}};
			l_ref_ = 1./160. * T_ * (pow(T_, 4) - 10. * pow(T_, 3) * sqrt(x0_[0]) + 80. * x0_[0] * (9. + x0_[0])
				+ 20. * pow(T_, 2) * (3. + 2. * x0_[0]) - 40. * T_ * sqrt(x0_[0]) * (9. + 2. * x0_[0]));
			g_ref_ = {-((T_ * (T_ - 4. * sqrt(x0_[0])) * (36. + pow(T_, 2) - 4. * T_ * sqrt(x0_[0]) + 8. * x0_[0])) / (32. * sqrt(x0_[0])))};
			H_ref_ = {{(T_ * (27. + pow(T_, 2) - 6. * T_ * sqrt(x0_[0]) + 12. * x0_[0])) / (12. * x0_[0])}};
		}


		template <
			typename VT1,
			typename VT2,
			typename VT3,
			typename VT4,
			typename VT5,
			typename MT1, bool SO1,
			typename MT2, bool SO2,
			typename MT3, bool SO3
		>
		void implicitDae(
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

			auto const z_hat = (~z)[0] + 1.;

			~f = {
				(~xdot)[0] + z_hat,
				(~x)[0] - pow(z_hat, 2),
			};
			
			~Jxdot = {
				{1.},
				{0.}
			};

			~Jx = {
				{0.},
				{1.}
			};

			~Jz = {
				{1.},
				{-2. * z_hat}
			};
		};


		template <
			typename VT1,
			typename VT2,
			typename MT1, bool SO1,
			typename VT3,
			typename VT4,
			typename MT2, bool SO2>
		void implicitDaeSensitivity(double t,
			blaze::Vector<VT1, blaze::columnVector> const& xdot,
			blaze::Vector<VT2, blaze::columnVector> const& x,
			blaze::Matrix<MT1, SO1> const& Sx,
			blaze::Vector<VT3, blaze::columnVector> const& z,
			blaze::Vector<VT4, blaze::columnVector> const& u,
			blaze::Matrix<MT2, SO2>& Sf) const
		{
			if (size(xdot) != NX || size(x) != NX || size(z) != NZ || size(u) != NU)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector size"));

			resize(Sf, NX + NZ, NX + NU);
			~Sf = {
				{0.},
				{(~Sx)(0, 0)}
			};
		};


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

			auto const z_hat = (~z)[0] + 1.;

			~r = {
				(~x)[0],
				3. * z_hat
			};

			~Sr = {
				{(~Sx)(0, 0)},
				{3. * (~Sz)(0, 0)}
			};
		}
	};


	TEST_F(IrkSimpleDaeTest, testBackwardEuler)
	{
		testIntegrate(BackwardEulerMethod {}, 0., 0.007);
	}


	TEST_F(IrkSimpleDaeTest, testGaussLegendre2)
	{
		testIntegrate(GaussLegendreMethod {2}, 0., 1e-7);
	}


	TEST_F(IrkSimpleDaeTest, testGaussLegendre3)
	{
		testIntegrate(GaussLegendreMethod {3}, 0., 1e-12);
	}


	TEST_F(IrkSimpleDaeTest, testGaussLegendre3Sensitivities)
	{
		testIntegrateWithSensitivities(GaussLegendreMethod {3}, 0., 1e-12);
	}


	TEST_F(IrkSimpleDaeTest, testGaussLegendre3LeastSquaresLagrangeTerm)
	{
		testIntegrateLeastSquaresLagrangeTerm(GaussLegendreMethod {3}, 0., 1e-12);
	}
}