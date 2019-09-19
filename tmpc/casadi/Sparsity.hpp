#pragma once

#include <tmpc/SizeT.hpp>

#include <casadi/mem.h>

#include <stdexcept>


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
			    throw std::invalid_argument("Invalid sparsity pattern passed to Sparsity ctor");
        }


        /// @brief Number of rows
        size_t rows() const
        {
            return sparsity_[0];
        }


        /// @brief Number of columns
        size_t columns() const
        {
            return sparsity_[1];
        }


        /// @brief Number of structural non-zeros
        size_t nnz() const
        {
            return colind()[columns()];
        }


        casadi_int const * colind() const
        {
            return sparsity_ + 2;
        }
        

		casadi_int const * rowind() const 
        {
            return colind() + columns() + 1;
        }


    private:
        casadi_int const * const sparsity_;
    };
}