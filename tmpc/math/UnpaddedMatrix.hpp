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
        UnpaddedMatrix(UnpaddedMatrix const&) = delete;
        UnpaddedMatrix(UnpaddedMatrix &&) = default;


        UnpaddedMatrix(size_t m, size_t n)
        :   Base {new Real[m * n], m, n}
        {
        }


        ~UnpaddedMatrix()
        {
            delete[] this->data();
        }
        

        using Base::operator=;
    };
}