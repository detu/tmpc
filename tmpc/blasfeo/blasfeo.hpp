#pragma once

#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include <tmpc/SizeT.hpp>


namespace tmpc :: blasfeo
{
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
    inline float const& element(blasfeo_smat const& m, int i, int j)
    {
        return BLASFEO_SMATEL(&m, i, j);
    }


    template <typename Real>
    class CustomMatrix
    {
    public:
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


    private:
        using Traits = detail::BlasfeoTraits<Real>;

        typename Traits::mat_t mat_;
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