#pragma once

#include "Types.hpp"
#include <tmpc/matrix/StorageOrder.hpp>
#include "MatrixAssign.hpp"
#include "Matrix.hpp"
#include "InitializerList.hpp"
#include "EigenBase.hpp"

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor 
{
    template <typename Type, size_t M, size_t N, StorageOrder SO>
    struct StaticMatrix;

    template <typename Type, size_t M, size_t N, StorageOrder SO>
    struct EigenBaseSelector<StaticMatrix<Type, M, N, SO>>
    {
        static auto constexpr storageOrder = 
            (M > 1 && N == 1) 
            ? Eigen::ColMajor 
            : (M == 1 && N > 1 
                ? Eigen::RowMajor 
                : (SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor));

        using type = Eigen::Matrix<Type, M, N, storageOrder>;
    };

    template <typename Type, size_t M, size_t N, StorageOrder SO = defaultStorageOrder>
    struct StaticMatrix
    :   Matrix<StaticMatrix<Type, M, N, SO>, SO>
    ,   EigenBase<StaticMatrix<Type, M, N, SO>>
    {
        using OurEigenBase = EigenBase<StaticMatrix<Type, M, N, SO>>;
        typedef Type ElementType;

        StaticMatrix()
        {            
        }

        StaticMatrix(Type const& rhs)
        {            
            this->setConstant(rhs);
        }

        StaticMatrix(initializer_list<initializer_list<Type>> list)
        {
            assign(*this, list);
        }

        template <typename T>
        StaticMatrix(Eigen::MatrixBase<T> const& rhs)
        :   OurEigenBase(rhs)
        {        
        }

        StaticMatrix& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }
    };

    /*
    template <typename Type, size_t M, size_t N, StorageOrder SO>
    struct EigenType<StaticMatrix<Type, M, N, SO>>
    {
        typedef typename StaticMatrix<Type, M, N, SO>::Base type;
    };
    */
}