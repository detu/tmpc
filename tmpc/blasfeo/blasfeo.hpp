#pragma once

#include <blasfeo_common.h>

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
        };
    }


    inline double& element(blasfeo_dmat& m, size_t i, size_t j)
    {
        
    }


    template <typename Real>
    class MatrixBase
    {
    private:
        using Traits = detail::BlasfeoTraits<Real>;

        typename Traits::mat_t mat_;
    };
}