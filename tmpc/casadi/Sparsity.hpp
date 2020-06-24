#pragma once

#include <tmpc/SizeT.hpp>

#include <casadi/mem.h>

#include <tmpc/Exception.hpp>
#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>


namespace tmpc :: casadi
{
	/// @brief Object-oriented view of Compressed Column Storage sparsity patterns.
    class Sparsity
    {
    public:
        Sparsity(casadi_int const * sparsity)
        :   rows_(sparsity[0])
        ,   columns_(sparsity[1])
        ,   colind_ {sparsity + 2}
        ,   rowind_ {colind_ + columns_ + 1}
        {
            if (colind_[0] != 0)
			    TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid sparsity pattern"));
        }


        /// @brief Number of rows
        size_t rows() const noexcept
        {
            return rows_;
        }


        /// @brief Number of columns
        size_t columns() const noexcept
        {
            return columns_;
        }


        /// @brief Number of structural non-zeros
        size_t nnz() const noexcept
        {
            return colind_[columns_];
        }


        /// @brief Copy data from a compressed column storage to a Blaze matrix.
        ///
        /// When this issue https://bitbucket.org/blaze-lib/blaze/issues/190/custom-matrix-vector-equivalent-for-sparse
        /// is resolved, we probably can use a wrapper for the compressed column storage (CCS) format
        /// profided by Blaze. Before that, we use functions which copy data to and from CSS arrays.
        template <typename MT, bool SO>
        inline void decompress(casadi_real const * data, blaze::Matrix<MT, SO>& m) const
        {
            if ((~m).rows() != rows_ || (~m).columns() != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Matrix size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    (~m)(rowind_[ind], j) = data[ind];
            }
        }


        template <typename VT>
        inline void decompress(casadi_real const * data, blaze::Vector<VT, blaze::columnVector>& v) const
        {
            if (size(v) != rows_ || 1 != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    (~v)[rowind_[ind]] = data[ind];
            }
        }


        template <typename VT>
        inline void decompress(casadi_real const * data, blaze::Vector<VT, blaze::rowVector>& v) const
        {
            if (1 != rows_ || size(v) != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    (~v)[j] = data[ind];
            }
        }


        /// @brief Copy data from a Blaze matrix to a compressed column storage.
        ///
        /// When this issue https://bitbucket.org/blaze-lib/blaze/issues/190/custom-matrix-vector-equivalent-for-sparse
        /// is resolved, we probably can use a wrapper for the compressed column storage (CCS) format
        /// profided by Blaze. Before that, we use functions which copy data to and from CSS arrays.
        template <typename MT, bool SO>
        inline void compress(blaze::Matrix<MT, SO> const& m, casadi_real * data) const
        {
            if ((~m).rows() != rows_ || (~m).columns() != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Matrix size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    data[ind] = (~m)(rowind_[ind], j);
            }
        }


        template <typename VT>
        inline void compress(blaze::Vector<VT, blaze::columnVector> const& v, casadi_real * data) const
        {
            if (size(v) != rows_ || 1 != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    data[ind] = (~v)[rowind_[ind]];
            }
        }


        template <typename VT>
        inline void compress(blaze::Vector<VT, blaze::rowVector> const& v, casadi_real * data) const
        {
            if (1 != rows_ || size(v) != columns_)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

            casadi_int ind = 0;
            for (size_t j = 0; j < columns_; ++j)
            {
                for (; ind < colind_[j + 1]; ++ind)
                    data[ind] = (~v)[j];
            }
        }


    private:
        size_t const rows_;
        size_t const columns_;
        casadi_int const * colind_;
        casadi_int const * rowind_;
    };
}