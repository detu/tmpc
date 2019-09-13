#pragma once

#include <tmpc/blasfeo/BlasfeoApi.hpp>


namespace tmpc :: blasfeo
{
    /// @brief BLASFEO matrix base class
    ///
    template <typename Derived>
    class Matrix
    {
    public:
        Derived& operator~() noexcept
        {
            return static_cast<Derived&>(*this);
        }


        Derived const& operator~() const noexcept
        {
            return static_cast<Derived const&>(*this);
        }


        decltype(auto) operator()(size_t i, size_t j) noexcept
        {
            return element(~*this, i, j);
        }


        decltype(auto) operator()(size_t i, size_t j) const noexcept
        {
            return element(~*this, i, j);
        }
        

    protected:
        /// \brief Protected constructor to prevent direct instantiation of Matrix objects.
        Matrix()
        {
        }
    };
}