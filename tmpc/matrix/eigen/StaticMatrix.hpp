#pragma once

#include "Types.hpp"
#include "StorageOrder.hpp"
#include "MatrixAssign.hpp"
#include "Matrix.hpp"
#include "EigenBase.hpp"

#include "Eigen.hpp"

namespace tmpc 
{
    template <typename Type, size_t M, size_t N, bool SO>
    struct StaticMatrix;

    template <typename Type, size_t M, size_t N, bool SO>
    struct EigenBaseSelector<StaticMatrix<Type, M, N, SO>>
    {
        static bool constexpr storageOrder = 
            (M > 1 && N == 1) 
            ? Eigen::ColMajor 
            : (SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor);

        using type = Eigen::Matrix<Type, M, N, storageOrder>;
    };

    template <typename Type, size_t M, size_t N, bool SO = defaultStorageOrder>
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
    template <typename Type, size_t M, size_t N, bool SO>
    struct EigenType<StaticMatrix<Type, M, N, SO>>
    {
        typedef typename StaticMatrix<Type, M, N, SO>::Base type;
    };
    */
}