#pragma once 

#include "StorageOrder.hpp"
//#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "Matrix.hpp"
#include "Rand.hpp"

#include "Eigen.hpp"

namespace tmpc 
{
    template <typename Type, bool SO>
    struct DynamicMatrix;

    template <typename Type, bool SO>
    struct EigenBaseSelector<DynamicMatrix<Type, SO>>
    {
        using type = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, 
            SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>;
    };

    template <typename Type, bool SO = defaultStorageOrder>
    struct DynamicMatrix
    :   Matrix<DynamicMatrix<Type, SO>, SO>
	,	EigenBase<DynamicMatrix<Type, SO>>
    {
        using Base = Matrix<DynamicMatrix<Type, SO>, SO>;
		using OurEigenBase = EigenBase<DynamicMatrix<Type, SO>>;
    	using ElementType = Type;

        DynamicMatrix()
        {
        }

        explicit DynamicMatrix(size_t M, size_t N)
        :   OurEigenBase(M, N)
        {        
        }

	    explicit DynamicMatrix(size_t M, size_t N, Type const& val)
	    :   OurEigenBase(M, N)
	    {
	        this->setConstant(val);
	    }

	    DynamicMatrix(initializer_list<initializer_list<Type>> list)
	    :   OurEigenBase(list.size(), determineColumns(list))
	    {
	        assign(*this, list);
	    }

        template <typename T>
        DynamicMatrix(Eigen::MatrixBase<T> const& rhs)
        :   OurEigenBase(rhs)
        {        
        }

        DynamicMatrix& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }
    };

    template <typename MT, bool SO>
    struct Rand<DynamicMatrix<MT, SO>>
    {
        decltype(auto) generate(size_t M, size_t N)
        {
            return DynamicMatrix<MT, SO>::Random(M, N);
        }
    };

}