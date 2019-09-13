#pragma once

#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include <tmpc/SizeT.hpp>

#include <memory>


namespace tmpc :: blasfeo
{
    /// \brief Alignment for BLASFEO data arrays
    inline static size_t constexpr alignment()
    {
        return 64;
    }


    namespace detail
    {
        template <typename Real>
        struct BlasfeoTraits;


        template <>
        struct BlasfeoTraits<double>
        {
            using mat_t = blasfeo_dmat;
            using vec_t = blasfeo_dvec;

            static auto constexpr create_mat = blasfeo_create_dmat;
        };


        /// \bried Aligned uninitialized array allocation
        template <typename T>
        T * alignedAlloc(size_t n)
        {
            return reinterpret_cast<T *>(std::aligned_alloc(alignment(), n * sizeof(T)));
        }


        /// \brief Deleter for aligned data arrays
        struct AlignedDeleter
        {
            template <typename T>
            void operator()(T * data)
            {
                std::free(data);
            }
        };
    }


    /// \brief Matrix element access
    inline double& element(blasfeo_dmat& m, size_t i, size_t j)
    {
        return BLASFEO_DMATEL(&m, i, j);
    }


    /// \brief Matrix element access
    inline float& element(blasfeo_smat& m, int i, int j)
    {
        return BLASFEO_SMATEL(&m, i, j);
    }


    /// \brief Const matrix element access
    inline double const& element(blasfeo_dmat const& m, int i, int j)
    {
        return BLASFEO_DMATEL(&m, i, j);
    }


    /// \brief Const matrix element access
    inline float const& element(blasfeo_smat const& m, int i, int j)
    {
        return BLASFEO_SMATEL(&m, i, j);
    }


    template <typename Real>
    class CustomMatrix
    {
    public:
        /// \brief Create a 0-by-0 matrix.
        CustomMatrix()
        {
            Traits::create_mat(0, 0, &mat_, nullptr);
        }

        
        CustomMatrix(Real * data, size_t m, size_t n)
        {
            Traits::create_mat(m, n, &mat_, data);
        }


        /// \brief Number of rows
        size_t rows() const
        {
            return mat_.m;
        }


        /// \brief Number of columns
        size_t columns() const
        {
            return mat_.n;
        }


        /// \brief Set new pointer and dimensions.
        /// Use with care!
        void reset(Real * data, size_t m, size_t n)
        {
            Traits::create_mat(m, n, &mat_, data);
        }


    private:
        using Traits = detail::BlasfeoTraits<Real>;

        typename Traits::mat_t mat_;
    };


    template <typename Real>
    class DynamicMatrix
    :   public CustomMatrix<Real>
    {
    public:
        DynamicMatrix(size_t m, size_t n)
        :   data_(detail::alignedAlloc<Real>(m * n))
        {
            this->reset(data_.get(), m, n);
        }


    private:
        std::unique_ptr<Real[], detail::AlignedDeleter> data_;
    };


    /// \brief Number of rows
    template <typename Real>
    size_t rows(CustomMatrix<Real> const& m)
    {
        return m.rows();
    }


    /// \brief Number of columns
    template <typename Real>
    size_t columns(CustomMatrix<Real> const& m)
    {
        return m.columns();
    }
}