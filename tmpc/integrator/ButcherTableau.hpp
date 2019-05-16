#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/SizeT.hpp>

#include <stdexcept>


namespace tmpc
{
	/**
	 * @brief Butcher tableau.
	 * @ingroup integrators
	 */
    template <typename Real>
    class ButcherTableau
    {
    public:
        using A_type = blaze::DynamicMatrix<Real>;
        using b_type = blaze::DynamicVector<Real, blaze::rowVector>;
        using c_type = blaze::DynamicVector<Real, blaze::columnVector>;


        ButcherTableau(A_type&& A, b_type&& b, c_type && c)
        :   A_ {std::move(A)}
        ,   b_ {std::move(b)}
        ,   c_ {std::move(c)}
        {
            if (!(rows(A_) == columns(A_) && rows(A_) == size(c_) && columns(A_) == size(b_)))
                throw std::invalid_argument(std::string(__func__) + ": inconsistent argumnent sizes");
        }


        A_type const& A() const
        {
            return A_;
        }


        b_type const& b() const
        {
            return b_;
        }


        c_type const& c() const
        {
            return c_;
        }


    private:
        A_type A_;
        b_type b_;
        c_type c_;
    };


    /// @brief Size of a Butcher tableau.
    template <typename Real>
    inline auto size(ButcherTableau<Real> const& t)
    {
        return size(t.b());
    }
}