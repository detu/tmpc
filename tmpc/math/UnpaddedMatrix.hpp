#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    template <typename Real, bool SO>
    class UnpaddedMatrix
    :   public blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>
    {
        using Base = blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>;


    public:
        UnpaddedMatrix()
        :   UnpaddedMatrix(0, 0)
        {
        }


        UnpaddedMatrix(UnpaddedMatrix const&) = delete;
        UnpaddedMatrix(UnpaddedMatrix &&) = default;


        UnpaddedMatrix(size_t m, size_t n)
        :   Base {new Real[m * n], m, n}
        {
        }


        template <typename MT, bool SO1>
        UnpaddedMatrix(blaze::Matrix<MT, SO1> const& rhs)
        :   UnpaddedMatrix(rows(rhs), columns(rhs))
        {
            *this = rhs;
        }


        ~UnpaddedMatrix()
        {
            delete[] this->data();
        }
        

        using Base::operator=;
    };
}