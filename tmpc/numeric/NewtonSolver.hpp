#pragma once

#include <tmpc/Math.hpp>

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
        ,   r_(nx)
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
            residualMaxNorm_ = inf<Real>();

			for (iterations_ = 0; iterations_ <= maxIterations_; ++iterations_)
			{
                fun(x_, r_, J_);
                residualMaxNorm_ = maxNorm(r_);

                monitor(iterations_, std::as_const(x_), std::as_const(r_), std::as_const(J_));

                // Residual within tolerance; exit the loop.
				if (residualMaxNorm_ < residualTolerance_)
					break;

                // Netwon method update: x(n+1) = x(n) - inv(J(n))*r(n)
                gesv(J_, r_, ipiv_.get());
                x_ -= r_;
            }

            if (!(residualMaxNorm_ < residualTolerance_))
                throw std::runtime_error("tmpc::NewtonSolver::solve(): max number of iteration reached but solution not found");

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
                throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": residual tolerance must be non-negative");

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


    private:
        size_t nx_;
        blaze::DynamicVector<Real, blaze::columnVector> x_;
        blaze::DynamicVector<Real, blaze::columnVector> r_;

        // The J matrix must be column-major in order blaze::gesv() to work.
        blaze::DynamicMatrix<Real, blaze::columnMajor> J_;

        size_t iterations_ = 0;
		size_t maxIterations_ = 10;
		Real residualMaxNorm_ = inf<Real>();
        Real residualTolerance_ = 1e-10;

        std::unique_ptr<int[]> ipiv_;
    };
}