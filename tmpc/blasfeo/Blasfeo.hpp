#pragma once

#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_s_aux.h>

#include <tmpc/SizeT.hpp>

#include <memory>


namespace tmpc :: blasfeo
{
    /// @brief BLASFEO matrix type selector
    template <typename Real>
    struct BlasfeoMatrix;


    template <>
    struct BlasfeoMatrix<double>
    {
        using type = blasfeo_dmat;
    };


    template <>
    struct BlasfeoMatrix<float>
    {
        using type = blasfeo_smat;
    };


    /// @brief BLASFEO vector type selector
    template <typename Real>
    struct BlasfeoVector;


    template <>
    struct BlasfeoVector<double>
    {
        using type = blasfeo_dvec;
    };


    template <>
    struct BlasfeoVector<float>
    {
        using type = blasfeo_svec;
    };


    /// \brief Number of rows
    size_t rows(blasfeo_dmat const& m)
    {
        return m.m;
    }


    /// \brief Number of rows
    size_t rows(blasfeo_smat const& m)
    {
        return m.m;
    }


    /// \brief Number of columns
    size_t columns(blasfeo_dmat const& m)
    {
        return m.n;
    }


    /// \brief Number of columns
    size_t columns(blasfeo_smat const& m)
    {
        return m.n;
    }


    /// \brief Alignment for BLASFEO data arrays
    inline static size_t constexpr alignment()
    {
        return 0x40;
    }


    // @brief returns the memory size (in bytes) needed for a BLASFEO matrix data with elements of type Real
    template <typename Real>
    size_t memsize_mat(size_t m, size_t n);


    // @brief returns the memory size (in bytes) needed for a double-precision BLASFEO matrix data
    template <>
    inline size_t memsize_mat<double>(size_t m, size_t n)
    {
        return blasfeo_memsize_dmat(m, n);
    }


    // @brief returns the memory size (in bytes) needed for a double-precision BLASFEO matrix data
    template <>
    inline size_t memsize_mat<float>(size_t m, size_t n)
    {
        return blasfeo_memsize_dmat(m, n);
    }


    /// @brief Create double-precision BLASFEO matrix
    inline void create_mat(size_t m, size_t n, blasfeo_dmat& sA, void * memory)
    {
        blasfeo_create_dmat(m, n, &sA, memory);
    }


    /// @brief Create single-precision BLASFEO matrix
    inline void create_mat(size_t m, size_t n, blasfeo_smat& sA, void * memory)
    {
        blasfeo_create_smat(m, n, &sA, memory);
    }


    namespace detail
    {
        /// \bried Aligned uninitialized array allocation
        template <typename T>
        T * alignedAlloc(size_t n)
        {
            return reinterpret_cast<T *>(std::aligned_alloc(alignment(), n * sizeof(T)));
        }


        /// \bried Aligned uninitialized array allocation
        template <>
        void * alignedAlloc<void>(size_t n)
        {
            return std::aligned_alloc(alignment(), n);
        }


        /// \brief Deleter for aligned data arrays
        struct AlignedDeleter
        {
            void operator()(void * data)
            {
                std::free(data);
            }
        };
    }


    /// @brief BLASFEO matrix type selector alias
    template <typename Real>
    using BlasfeoMatrix_t = typename BlasfeoMatrix<Real>::type;


    /// @brief BLASFEO vector type selector alias
    template <typename Real>
    using BlasfeoVector_t = typename BlasfeoVector<Real>::type;


    /// \brief Matrix element access
    inline double& element(blasfeo_dmat& m, size_t i, size_t j)
    {
        return BLASFEO_DMATEL(&m, i, j);
    }


    /// \brief Matrix element access
    inline float& element(blasfeo_smat& m, size_t i, size_t j)
    {
        return BLASFEO_SMATEL(&m, i, j);
    }


    /// \brief Const matrix element access
    inline double const& element(blasfeo_dmat const& m, size_t i, size_t j)
    {
        return BLASFEO_DMATEL(&m, i, j);
    }


    /// \brief Const matrix element access
    inline float const& element(blasfeo_smat const& m, size_t i, size_t j)
    {
        return BLASFEO_SMATEL(&m, i, j);
    }


    /// @brief BLASFEO matrix base class
    ///
    template <typename Real>
    class MatrixBase
    :   public BlasfeoMatrix_t<Real>
    {
    protected:
        /// \brief Protected constructor to prevent direct instantiation of MatrixBase objects.
        MatrixBase()
        {
        }
    };



    /// @brief BLASFEO matrix using preallocated memory array
    ///
    template <typename Real>
    class CustomMatrix
    :   public MatrixBase<Real>
    {
    public:
        /// \brief Create a 0-by-0 matrix.
        CustomMatrix()
        {
            create_mat(0, 0, this, nullptr);
        }

        
        CustomMatrix(Real * data, size_t m, size_t n)
        {
            create_mat(m, n, *this, data);
        }


        /// \brief Set new pointer and dimensions.
        /// Use with care!
        void reset(Real * data, size_t m, size_t n)
        {
            create_mat(m, n, *this, data);
        }
    };


    /// @brief BLASFEO matrix managing its own memory pool
    ///
    template <typename Real>
    class DynamicMatrix
    :   public MatrixBase<Real>
    {
    public:
        /// \brief Create a 0-by-0 matrix.
        DynamicMatrix()
        {
            create_mat(0, 0, *this, nullptr);
        }


        DynamicMatrix(size_t m, size_t n)
        :   data_(reinterpret_cast<Real *>(std::aligned_alloc(alignment(), memsize_mat<Real>(m, n))))
        {
            create_mat(m, n, *this, data_.get());
        }


    private:
        std::unique_ptr<Real[], detail::AlignedDeleter> data_;
    };
}