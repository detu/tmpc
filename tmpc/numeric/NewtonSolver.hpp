#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/Exception.hpp>

#include <memory>


namespace tmpc
{
    /// @brief Newton method solver for systems of non-linear equations.
    template <typename Real>
    class NewtonSolver
    {
    public:
        explicit NewtonSolver(size_t nx)
        :   nx_(nx)
        ,   x1_(nx)
        ,   r_(nx)
        ,   r1_(nx)
        ,   d_(nx)
        ,   J_(nx, nx, Real {})
        ,   ipiv_(new int[nx])
        {
        }


        template <typename F, typename VT1, typename VT2>
        void operator()(F const& fun,
            blaze::Vector<VT1, blaze::columnVector> const& x0,
            blaze::Vector<VT2, blaze::columnVector>& xf)
        {
            (*this)(fun, x0, xf, emptyMonitor());
        }


        template <typename F, typename VT1, typename VT2, typename Monitor>
        void operator()(F const& fun,
            blaze::Vector<VT1, blaze::columnVector> const& x0,
            blaze::Vector<VT2, blaze::columnVector>& xf,
            Monitor&& monitor)
        {
            if (size(x0) != nx_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector dimension"));

            functionEvaluations_ = 0;

            ~xf = ~x0;
            fun(~xf, r_, J_);
            factorizeJacobian();
            ++functionEvaluations_;

            for (iterations_ = 0; iterations_ <= maxIterations_; ++iterations_)
			{
                residualMaxNorm_ = maxNorm(r_);

                monitor(iterations_, std::as_const(~xf), std::as_const(r_), std::as_const(J_));

                // Residual within tolerance; exit the loop.
				if (residualMaxNorm_ < residualTolerance_)
					break;

                // Calculate search direction d(n)=-inv(J(n))*r(n)
                jacobianSolve(r_, d_);

                // Step size
                Real t = 1.;

                // Do backtracking search.
                // alpha == 1. disables backtracking.
                while (fun(x1_ = ~xf + t * d_, r1_, J_), ++functionEvaluations_, alpha_ < 1. && !allAbsLessThan(r1_, r_))
                    t *= alpha_;
                factorizeJacobian();
                
                // Netwon method update: x(n+1) = x(n) + t*d
                ~xf = x1_;
                r_ = r1_;
            }

            if (!(residualMaxNorm_ < residualTolerance_))
                TMPC_THROW_EXCEPTION(std::runtime_error("Max number of iteration reached but solution not found"));
        }


        /// @brief Solve the equation using Newton method and calculate 
        /// solution sensitivities w.r.t. parameters.
        ///
        /// @param dxf_dp a NX-by-NP matrix of solution sensitivities
        ///
        template <
            typename F, 
            typename DFDP, 
            typename VT1, 
            typename VT2, 
            typename MT
        >
        void operator()(F const& fun,
            DFDP const& dfdp,
            blaze::Vector<VT1, blaze::columnVector> const& x0,
            blaze::DenseVector<VT2, blaze::columnVector>& xf,
            blaze::DenseMatrix<MT, blaze::columnMajor>& dxf_dp)
        {
            if (size(x0) != nx_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid vector dimension"));
                
            (*this)(fun, ~x0, ~xf);

            // Calculate df(x,p)/dp at the solution
            dfdp(~xf, ~dxf_dp);

            // From 0 = df(x,p)/dx * dx^*/dp + df(x,p)/dp
            // calculate dx^*/dp = -inv(df(x,p)/dx) * df(x,p)/dp
            jacobianSolve(~dxf_dp, ~dxf_dp);
        }


        size_t maxIterations() const noexcept
        {
            return maxIterations_;
        }


        void maxIterations(size_t val)
        {
            maxIterations_ = val;
        }


        Real residualTolerance() const noexcept
        {
            return residualTolerance_;
        }


        void residualTolerance(Real val)
        {
            if (val < 0)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Residual tolerance must be non-negative"));

            residualTolerance_ = val;
        }


        Real residualMaxNorm() const noexcept
        {
            return residualMaxNorm_;
        }


        /// @brief Total number of Newton iterations during last solve.
        size_t iterations() const noexcept
        {
            return iterations_;
        }


        /// @brief Total number of function evaluations during last solve.
        size_t functionEvaluations() const noexcept
        {
            return functionEvaluations_;
        }


        Real backtrackingAlpha() const noexcept
        {
            return alpha_;
        }


        /// @brief Set backtracking alpha parameter value.
        ///
        /// alpha must be withing the range (0., 1.].
        /// Setting alpha = 1. disables backtracking.
        void backtrackingAlpha(Real val)
        {
            if (!(0. < val && val <= 1.))
                TMPC_THROW_EXCEPTION(std::invalid_argument("Backtracking alpha must be within the (0, 1] range"));

            alpha_ = val;
        }


        /// @brief Solve J*x+b=0 w.r.t. x where x, b are vectors and J is the current Jacobian.
        template <typename VT1, typename VT2>
        void jacobianSolve(
            blaze::DenseVector<VT1, blaze::columnVector> const& b, 
            blaze::DenseVector<VT2, blaze::columnVector>& x) const
        {
            ~x = -~b;
            getrs(~J_, ~x, 'N', ipiv_.get());
        }


        /// @brief Solve J*X+B=0 w.r.t. x where X, B are matrices and J is the current Jacobian.
        template <typename MT1, typename MT2>
        void jacobianSolve(
            blaze::DenseMatrix<MT1, blaze::columnMajor> const& B, 
            blaze::DenseMatrix<MT2, blaze::columnMajor>& X) const
        {
            ~X = -~B;
            getrs(~J_, ~X, 'N', ipiv_.get());
        }


    private:
        size_t nx_;
        blaze::DynamicVector<Real, blaze::columnVector> x1_;
        blaze::DynamicVector<Real, blaze::columnVector> r_;
        blaze::DynamicVector<Real, blaze::columnVector> r1_;
        blaze::DynamicVector<Real, blaze::columnVector> d_;

        // Factorized Jacobian.
        // Must be column-major in order blaze::getrf() and blaze::getrs() to work as expected.
        blaze::DynamicMatrix<Real, blaze::columnMajor> J_;

        size_t iterations_ = 0;
		size_t maxIterations_ = 10;
        size_t functionEvaluations_ = 0;
		Real residualMaxNorm_ = inf<Real>();
        Real residualTolerance_ = 1e-10;

        // Backtracking alpha
        Real alpha_ = 1.;

        std::unique_ptr<int[]> ipiv_;


        template <typename VT1, typename VT2, bool TF>
        static bool allAbsLessThan(blaze::Vector<VT1, TF> const& a, blaze::Vector<VT2, TF> const& b)
        {
            size_t const n = size(a);

            if (n != size(b))
                TMPC_THROW_EXCEPTION(std::invalid_argument("Vector sizes don't match"));

            size_t i = 0; 
            while (i < n && abs((~a)[i]) < abs((~b)[i]))
                ++i;
                
            return i >= n;
        }


        void factorizeJacobian()
        {
            // Factorize the Jacobian
            getrf(J_, ipiv_.get());
        }


        static auto emptyMonitor()
        {
            return [] (size_t, auto const&, auto const&, auto const&) {};
        }
    };
}