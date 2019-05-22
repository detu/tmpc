#pragma once

#include <tmpc/Math.hpp>

#include <boost/throw_exception.hpp>

#include <memory>
#include <stdexcept>


namespace tmpc
{
    /// @brief Newton method solver for systems of non-linear equations.
    template <typename Real>
    class NewtonSolver
    {
    public:
        NewtonSolver(size_t nx)
        :   nx_(nx)
        ,   x_(nx)
        ,   x1_(nx)
        ,   r_(nx)
        ,   r1_(nx)
        ,   d_(nx)
        ,   J_(nx, nx, Real {})
        ,   ipiv_(new int[nx])
        {
        }


        template <typename F, typename VT>
        auto const& solve(F const& fun, blaze::Vector<VT, blaze::columnVector> const& x0)
        {
            return solve(fun, x0, [] (size_t, auto const&, auto const&, auto const&) {});
        }


        template <typename F, typename VT, typename Monitor>
        auto const& solve(F const& fun, blaze::Vector<VT, blaze::columnVector> const& x0, Monitor monitor)
        {
            x_ = x0;
            fun(x_, r_, J_);

			for (iterations_ = 0; iterations_ <= maxIterations_; ++iterations_)
			{
                residualMaxNorm_ = maxNorm(r_);

                monitor(iterations_, std::as_const(x_), std::as_const(r_), std::as_const(J_));

                // Residual within tolerance; exit the loop.
				if (residualMaxNorm_ < residualTolerance_)
					break;

                // Calculate search direction d(n)=-inv(J(n))*r(n)
                d_ = -r_;
                gesv(J_, d_, ipiv_.get());

                // Step size
                Real t = 1.;

                // Do backtracking search
                while (fun(x1_ = x_ + t * d_, r1_, J_), !allAbsLessThan(r1_, r_) /*maxNorm(r1_) >= residualMaxNorm_*/)
                    t *= alpha_;

                // Netwon method update: x(n+1) = x(n) + t*d
                x_ = x1_;
                r_ = r1_;
            }

            if (!(residualMaxNorm_ < residualTolerance_))
                BOOST_THROW_EXCEPTION(std::runtime_error("Max number of iteration reached but solution not found"));

            return x_;
        }


        size_t maxIterations() const
        {
            return maxIterations_;
        }


        void maxIterations(size_t val)
        {
            maxIterations_ = val;
        }


        Real residualTolerance() const
        {
            return residualTolerance_;
        }


        void residualTolerance(Real val)
        {
            if (val < 0)
                BOOST_THROW_EXCEPTION(std::invalid_argument("Residual tolerance must be non-negative"));

            residualTolerance_ = val;
        }


        Real residualMaxNorm() const
        {
            return residualMaxNorm_;
        }


        size_t iterations() const
        {
            return iterations_;
        }


        Real backTrackingAlpha() const
        {
            return alpha_;
        }


        void backTrackingAlpha(Real val)
        {
            if (!(0. < val && val < 1.))
                throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": backtracking alpha must be within (0, 1)");

            alpha_ = val;
        }


    private:
        size_t nx_;
        blaze::DynamicVector<Real, blaze::columnVector> x_;
        blaze::DynamicVector<Real, blaze::columnVector> x1_;
        blaze::DynamicVector<Real, blaze::columnVector> r_;
        blaze::DynamicVector<Real, blaze::columnVector> r1_;
        blaze::DynamicVector<Real, blaze::columnVector> d_;

        // The J matrix must be column-major in order blaze::gesv() to work.
        blaze::DynamicMatrix<Real, blaze::columnMajor> J_;

        size_t iterations_ = 0;
		size_t maxIterations_ = 10;
		Real residualMaxNorm_ = inf<Real>();
        Real residualTolerance_ = 1e-10;

        // Backtracking alpha
        Real alpha_ = 0.5;

        std::unique_ptr<int[]> ipiv_;


        template <typename VT1, typename VT2, bool TF>
        static bool allAbsLessThan(blaze::Vector<VT1, TF> const& a, blaze::Vector<VT2, TF> const& b)
        {
            size_t const n = size(a);

            if (n != size(b))
                throw std::invalid_argument(std::string(__func__) + ": vector sizes don't match");

            size_t i = 0; 
            while (i < n && abs((~a)[i]) < abs((~b)[i]))
                ++i;
                
            return i >= n;
        }
    };
}