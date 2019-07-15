#pragma once

#include <tmpc/Math.hpp>

#include <boost/throw_exception.hpp>

#include <iosfwd>
#include <stdexcept>


namespace tmpc
{
    /// @brief Discrete-time Algebraic Riccati equation (DARE) solver.
    template <typename Real>
    class Dare
    {
    public:
        /// @brief Constructor
        ///
        /// @param nx number of states
        /// @param nu number of controls
        Dare(size_t nx, size_t nu)
        :   nx_(nx)
        ,   nu_(nu)
        ,   H_(2 * nx + nu, 2 * nx + nu, Real {})
        ,   J_(2 * nx + nu, 2 * nx + nu, Real {})
        ,   QR_q_(2 * nx + nu, 2 * nx + nu)
        ,   QR_r_(nu, nu)
        ,   H_compr_(2 * nx, 2 * nx)
        ,   J_compr_(2 * nx, 2 * nx)
        ,   QZ_q_(2 * nx, 2 * nx)
        ,   QZ_z_(2 * nx, 2 * nx)
        ,   QZ_alpha_(2 * nx)
        ,   QZ_beta_(2 * nx)
        ,   ipiv_(new int[nx])
        {
        }


        /// @brief Solves DARE of the form A^T*X*A−X−A^T*X*B*(B^T*X*B+R)^{−1}B^T*X*A+Q=0
        ///
        /// TODO: should Q, R, X be declared as a SymmetricMatrix<> ?
        template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3, typename MT4, bool SO4,
            typename MT5, bool SO5>
        void operator()(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& B,
            blaze::Matrix<MT3, SO3> const& Q, blaze::Matrix<MT4, SO4> const& R,
            blaze::Matrix<MT5, SO5>& X)
        {
            //                    [  A   0   B  ]       [ E   0   0 ]
            //        H - t J  =  [ -Q   E' -S  ]  - t  [ 0   A'  0 ]
            //                    [  S'  0   R  ]       [ 0  -B'  0 ]

            submatrix(H_, 0, 0, nx_, nx_) = ~A;
            submatrix(H_, 0, 2 * nx_, nx_, nu_) = ~B;
            submatrix(H_, nx_, 0, nx_, nx_) = -(~Q);
            submatrix(H_, nx_, nx_, nx_, nx_) = blaze::IdentityMatrix<Real>(nx_);
            // submatrix(H_, nx_, 2 * nx_, nx_, nu_) = -(~S);
            // submatrix(H_, 2 * nx_, 0, nu_, nx_) = trans(~S);
            submatrix(H_, 2 * nx_, 2 * nx_, nu_, nu_) = ~R;

            submatrix(J_, 0, 0, nx_, nx_) = blaze::IdentityMatrix<Real>(nx_);
            submatrix(J_, nx_, nx_, nx_, nx_) = trans(~A);
            submatrix(J_, 2 * nx_, nx_, nu_, nx_) = -trans(~B);

            // Compression step on H(:,n2+1:n2+m) = [S1;-S2;R]
            if (nu_ > 0)
            {
                // [q,r] = qr(H(:,n2+1:n2+m));
                qr(submatrix(H_, 0, 2 * nx_, rows(H_), nu_), QR_q_, QR_r_, 2 * nx_ + nu_);
                // H = q(:,m+1:n2+m)'*H(:,1:n2);
                H_compr_ = trans(submatrix(QR_q_, 0, nu_, rows(QR_q_), 2 * nx_)) * submatrix(H_, 0, 0, rows(H_), 2 * nx_);
                // J = q(:,m+1:n2+m)'*J(:,1:n2);
                J_compr_ = trans(submatrix(QR_q_, 0, nu_, rows(QR_q_), 2 * nx_)) * submatrix(J_, 0, 0, rows(J_), 2 * nx_);
            }

            // Perform Schur decomposition and reorder matrix elements and eigenvalues
            // such that the eigeivalues outside the unit disk go first.
            gges(&outsideUnitDisk, J_compr_, H_compr_, QZ_alpha_, QZ_beta_, QZ_q_, QZ_z_);

            auto X1 = submatrix(QZ_z_, 0, 0, nx_, nx_);
            auto X2 = submatrix(QZ_z_, nx_, 0, nx_, nx_);

            areCheckout(X1, X2, QZ_alpha_, QZ_beta_, X);

            // Solve X * X1 = X2
            X1 = trans(X1);
            X2 = trans(X2);
            gesv(X1, X2, ipiv_.get());

            // Symmetrize
            ~X = 0.5 * (X2 + trans(X2));
        }


        void dump(std::ostream& os) const
        {
            os << "H_ = \n" << H_ << "\n";
            os << "J_ = \n" << J_ << "\n";
            os << "QR_q_ = \n" << QR_q_ << "\n";
            os << "QR_r_ = \n" << QR_r_ << "\n";
            os << "H_compr_ = \n" << H_compr_ << "\n";
            os << "J_compr_ = \n" << J_compr_ << "\n";
            os << "QZ_q_ = \n" << QZ_q_ << "\n";
            os << "QZ_z_ = \n" << QZ_z_ << "\n";
            os << "abs(lambda) = " << abs(QZ_alpha_ / QZ_beta_) << "\n";
        }


    private:
        size_t const nx_;
        size_t const nu_;
        blaze::DynamicMatrix<Real> H_;
        blaze::DynamicMatrix<Real> J_;
        blaze::DynamicMatrix<Real> QR_q_;
        blaze::DynamicMatrix<Real> QR_r_;
        blaze::DynamicMatrix<Real, blaze::columnMajor> H_compr_;
        blaze::DynamicMatrix<Real, blaze::columnMajor> J_compr_;
        blaze::DynamicMatrix<Real, blaze::rowMajor> QZ_q_;
        blaze::DynamicMatrix<Real, blaze::columnMajor> QZ_z_;
        blaze::DynamicVector<std::complex<Real>, blaze::rowVector> QZ_alpha_;
        blaze::DynamicVector<Real, blaze::rowVector> QZ_beta_;
        std::unique_ptr<int> ipiv_;


        static int outsideUnitDisk(Real const * alphar, Real const * alphai, Real const * beta)
        {
            return outsideUnitDisk({*alphar, *alphai}, *beta);
        }


        static bool outsideUnitDisk(std::complex<Real> alpha, Real beta)
        {
            return pow(real(alpha), 2) + pow(imag(alpha), 2) > pow(beta, 2);
        }


        /// @brief Checks for proper extraction of stable invariant subspace.
        template <typename MT1, bool SO1, typename MT2, bool SO2, 
            typename VT1, bool TF1, typename VT2, bool TF2, typename MT3, bool SO3>
        void areCheckout(blaze::Matrix<MT1, SO1> const& X1, blaze::Matrix<MT2, SO2> const& X2,
            blaze::Vector<VT1, TF1> const& alpha, blaze::Vector<VT2, TF2> const& beta,
            blaze::Matrix<MT3, SO3>& X12) const
        {
            ~X12 = trans(~X1) * (~X2);
            auto const asym = l1Norm(~X12 - trans(~X12));

            // Check that first nx_ eigenvalues are outside the unit disk, and the other nx_ are inside the unit disk.
            bool eigenvalues_ok = true;
            for (size_t i = 0; i < 2 * nx_ && eigenvalues_ok; ++i)
                eigenvalues_ok = eigenvalues_ok && (outsideUnitDisk((~alpha)[i], (~beta)[i]) == i < nx_);

            // Check solution asymmetry
            if (!eigenvalues_ok || asym > std::max(1e-3 * std::numeric_limits<Real>::epsilon(), 0.1 * l1Norm(~X12)))
                BOOST_THROW_EXCEPTION(std::runtime_error(
                    "Could not (reliably) isolate stable invariant subspace of dimension " + std::to_string(nx_)));

            // Check accuracy
            if (asym > sqrt(std::numeric_limits<Real>::epsilon()))
                BOOST_THROW_EXCEPTION(std::runtime_error("Riccati equation solution accuracy is not sufficient"));
        }
    };
}