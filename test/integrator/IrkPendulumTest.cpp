#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>

#include <casadi/casadi.hpp>

#include <stdexcept>


namespace tmpc :: testing
{
	using Real = double;


	template <typename VT, bool TF>
	inline casadi::DM toDM(blaze::DenseVector<VT, TF> const& v)
	{
		casadi::DM dm(casadi::Sparsity::dense(
			TF == blaze::columnVector ? size(v) : 1, TF == blaze::rowVector ? size(v) : 1));
		std::copy_n(data(v), size(v), dm.ptr());

		return dm;
	}


	template <bool TF>
	inline auto asVector(casadi::DM const& dm)
	{
		if (TF == blaze::columnVector)
		{
			if (dm.size2() != 1)
				throw std::invalid_argument("casadi::DM value is not a column vector");
		}
		else
		{
			if (dm.size1() != 1)
				throw std::invalid_argument("casadi::DM value is not a row vector");
		}

		return blaze::CustomVector<double const, blaze::unaligned, blaze::unpadded, TF>(
			dm.ptr(), TF == blaze::columnVector ? dm.size1() : dm.size2());
	}


	inline auto asMatrix(casadi::DM const& dm)
	{
		if (!dm.is_dense())
			throw std::invalid_argument("Sparse casadi::DM are not supported");

		return blaze::CustomMatrix<double const, blaze::unaligned, blaze::unpadded, blaze::columnMajor>(
			dm.ptr(), dm.size1(), dm.size2());
	}

	
	class IrkPendulumTest 
	: 	public Test
	{
	protected:
		void SetUp() override
		{
			casadi::MX t = casadi::MX::sym("t");	// time
    		casadi::MX phi = casadi::MX::sym("phi");	// the angle phi
    		casadi::MX omega = casadi::MX::sym("omega");	// the first derivative of phi w.r.t time
    		casadi::MX u = casadi::MX::sym("u");	// force acting on the pendulum
    		double l = 1.;	// the length of the pendulum
    		double m = 1.;	// the mass of the pendulum
    		double g = 9.81;	// the gravitational constant
    		double alpha = 2.;	// friction constant
    		casadi::MX x = vertcat(phi, omega);	// state vector

			casadi::MX z = sin(phi);
    		casadi::MX dx = vertcat(omega, -(m*g/l)*z - alpha*omega + u/(m*l));
    		casadi::MX q = vertcat(pow(omega, 2), phi);	// quadrature term

			// Create integrator function
			integrator_ = casadi::integrator("pendulum_integrator", "cvodes", 
				casadi::MXDict {{"x", x}, {"p", u}, {"ode", dx}, {"quad", q}}, casadi::Dict {{"tf", h_}});


			// Create ODE function
			ode_ = casadi::Function("pendulum_ode", 
				casadi::MXVector {t, x, u}, 
                casadi::MXVector {densify(dx), densify(jacobian(dx, x)), densify(jacobian(dx, u)), 
                    densify(q), densify(jacobian(q, x)), densify(jacobian(q, u))},
                casadi::StringVector {"t", "x", "u"},
				casadi::StringVector {"xdot", "A", "B", "q", "qA", "qB"});
		}


		/// @brief Integrate using the reference integrator.
		blaze::DynamicVector<Real, blaze::columnVector> refIntegrate(
			double t0,
			blaze::DynamicVector<Real, blaze::columnVector> const& x0,
			blaze::DynamicVector<Real, blaze::columnVector> const& u
		)
		{
			if (t0 != 0.)
				throw std::invalid_argument("t0 other than 0 is currently not supported");

			casadi::DMDict out = integrator_(casadi::DMDict {{"x0", toDM(x0)}, {"p", toDM(u)}});
			return asVector<blaze::columnVector>(out["xf"]);
		}


		/// @brief Evaluate ODE.
		template <typename VT1, typename VT2, typename VT3, typename MT1, bool SO1, typename MT2, bool SO2>
		void ode(
			double t,
			blaze::DenseVector<VT1, blaze::columnVector> const& x,
			blaze::DenseVector<VT2, blaze::columnVector> const& u,
			blaze::Vector<VT3, blaze::columnVector>& xdot,
			blaze::Matrix<MT1, SO1>& A,
			blaze::Matrix<MT2, SO2>& B)
		{
			casadi::DMDict in {{"t", t}, {"x", toDM(x)}, {"u", toDM(u)}};
			casadi::DMDict out = ode_(in);

			~xdot = asVector<blaze::columnVector>(out.at("xdot"));
			~A = asMatrix(out.at("A"));
			~B = asMatrix(out.at("B"));
		}


		size_t nx() const
		{
			return 2;
		}


		size_t nz() const
		{
			return 0;
		}


		size_t nu() const
		{
			return 1;
		}


		template <typename Method>
		void testIntegrate(Method const& method, double abs_tol, double rel_tol)
		{
			ImplicitRungeKutta<double> irk(method, nx(), nz(), nu());

			double t0 = 0.;

			size_t const num_points = 1000;
			blaze::DynamicVector<double, blaze::columnVector> x0(nx());
			blaze::DynamicVector<double, blaze::columnVector> u(nu());

			for (size_t i = 0; i < num_points; ++i)
			{
				randomize(x0);
				randomize(u);

				blaze::DynamicVector<double> xf(nx());
				irk(
					[this] (double t, auto const& xdot, auto const& x, auto const& z,
						auto const& u, auto& f, auto& Jxdot, auto& Jx, auto& Jz)
					{
						blaze::DynamicMatrix<Real> df_du(nx(), nu());
						ode(t, x, u, f, Jx, df_du);

						f -= xdot;
						Jxdot = -blaze::IdentityMatrix<double>(nx());
					}, 
					t0, h_, x0, u, xf
				);

				blaze::DynamicVector<double, blaze::columnVector> const xf_ref = refIntegrate(t0, x0, u);

				TMPC_EXPECT_APPROX_EQ(xf, xf_ref, abs_tol, rel_tol);
			}
		}


	private:
		casadi::Function integrator_;
		casadi::Function ode_;
		double const h_ = 0.01;	// time step
	};


	TEST_F(IrkPendulumTest, testBackwardEuler)
	{
		testIntegrate(BackwardEulerMethod(), 1e-3, 2e-3);
	}


	TEST_F(IrkPendulumTest, testGaussLegendre2)
	{
		testIntegrate(GaussLegendreMethod(2), 1e-6, 1e-4);
	}


	TEST_F(IrkPendulumTest, testGaussLegendre3)
	{
		testIntegrate(GaussLegendreMethod(3), 1e-6, 1e-5);
	}
}