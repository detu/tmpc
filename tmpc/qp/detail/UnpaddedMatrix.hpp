#pragma once

#include <tmpc/Matrix.hpp>


namespace tmpc :: detail
{
    template <typename Kernel, StorageOrder SO>
    class UnpaddedMatrix
    :   public CustomMatrix<Kernel, unaligned, unpadded, SO>
    {
        using Base = CustomMatrix<Kernel, unaligned, unpadded, SO>;
        using Real = typename Kernel::Real;

    public:
        UnpaddedMatrix(size_t m, size_t n)
        :   CustomMatrix<Kernel, unaligned, unpadded, SO> {new Real[m * n], m, n}
        {
        }

        ~UnpaddedMatrix()
        {
            delete[] this->data();
        }

        using Base::operator=;
    };
}