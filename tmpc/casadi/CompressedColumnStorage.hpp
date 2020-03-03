#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/casadi/Sparsity.hpp>
#include <tmpc/Exception.hpp>

#include <casadi/mem.h>

#include <blaze/Math.h>


namespace tmpc :: casadi
{
	/// @brief Copy data from a compressed column storage to a Blaze matrix.
	///
	/// When this issue https://bitbucket.org/blaze-lib/blaze/issues/190/custom-matrix-vector-equivalent-for-sparse
	/// is resolved, we probably can use a wrapper for the compressed column storage (CCS) format
	/// profided by Blaze. Before that, we use functions which copy data to and from CSS arrays.
	template <typename MT, bool SO>
	inline void fromCompressedColumnStorage(casadi_real const * data, Sparsity const& sparsity, blaze::Matrix<MT, SO>& m)
	{
        if (rows(m) != sparsity.rows() || columns(m) != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Matrix size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				(~m)(sparsity.rowind()[ind], j) = data[ind];
		}
	}


	template <typename VT>
	inline void fromCompressedColumnStorage(casadi_real const * data, Sparsity const& sparsity, blaze::Vector<VT, blaze::columnVector>& v)
	{
        if (size(v) != sparsity.rows() || 1 != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				(~v)[sparsity.rowind()[ind]] = data[ind];
		}
	}


	template <typename VT>
	inline void fromCompressedColumnStorage(casadi_real const * data, Sparsity const& sparsity, blaze::Vector<VT, blaze::rowVector>& v)
	{
        if (1 != sparsity.rows() || size(v) != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				(~v)[j] = data[ind];
		}
	}


	/// @brief Copy data from a Blaze matrix to a compressed column storage.
	///
	/// When this issue https://bitbucket.org/blaze-lib/blaze/issues/190/custom-matrix-vector-equivalent-for-sparse
	/// is resolved, we probably can use a wrapper for the compressed column storage (CCS) format
	/// profided by Blaze. Before that, we use functions which copy data to and from CSS arrays.
	template <typename MT, bool SO>
	inline void toCompressedColumnStorage(blaze::Matrix<MT, SO> const& m, casadi_real * data, Sparsity const& sparsity)
	{
		if (rows(m) != sparsity.rows() || columns(m) != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Matrix size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				data[ind] = (~m)(sparsity.rowind()[ind], j);
		}
	}


	template <typename VT>
	inline void toCompressedColumnStorage(blaze::Vector<VT, blaze::columnVector> const& v, casadi_real * data, Sparsity const& sparsity)
	{
        if (size(v) != sparsity.rows() || 1 != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				data[ind] = (~v)[sparsity.rowind()[ind]];
		}
	}


	template <typename VT>
	inline void toCompressedColumnStorage(blaze::Vector<VT, blaze::rowVector> const& v, casadi_real * data, Sparsity const& sparsity)
	{
        if (1 != sparsity.rows() || size(v) != sparsity.columns())
			TMPC_THROW_EXCEPTION(std::invalid_argument("Vector size does not match the sparsity pattern"));

		casadi_int ind = 0;
		for (size_t j = 0; j < sparsity.columns(); ++j)
		{
			for (; ind < sparsity.colind()[j + 1]; ++ind)
				data[ind] = (~v)[j];
		}
	}
}