#pragma once

#include <tmpc/SizeT.hpp>

#include <casadi/mem.h>

#include <tmpc/Exception.hpp>


namespace tmpc :: casadi
{
	/// @brief Object-oriented view of Compressed Column Storage sparsity patterns.
    class Sparsity
    {
    public:
        Sparsity(casadi_int const * sparsity)
        :   sparsity_(sparsity)
        {
            if (colind()[0] != 0)
			    TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid sparsity pattern"));
        }


        /// @brief Number of rows
        size_t rows() const noexcept
        {
            return sparsity_[0];
        }


        /// @brief Number of columns
        size_t columns() const noexcept
        {
            return sparsity_[1];
        }


        /// @brief Number of structural non-zeros
        size_t nnz() const noexcept
        {
            return colind()[columns()];
        }


        casadi_int const * colind() const noexcept
        {
            return sparsity_ + 2;
        }
        

		casadi_int const * rowind() const noexcept
        {
            return colind() + columns() + 1;
        }


    private:
        casadi_int const * const sparsity_;
    };
}