#pragma once 

#include "StorageOrder.hpp"
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "Matrix.hpp"

#include <Eigen/Dense>

namespace tmpc 
{
    template <typename Type, bool SO = defaultStorageOrder>
    struct DynamicMatrix
    :   Matrix<DynamicMatrix<Type, SO>, SO>
	,	EigenBase<DynamicMatrix<Type, SO>>
    {
        using Base = Matrix<DynamicMatrix<Type, SO>, SO>;
		using OurEigenBase = EigenBase<DynamicMatrix<Type, SO>>;
    	using ElementType = Type;

        DynamicMatrix(size_t M, size_t N)
        :   OurEigenBase(M, N)
        {        
        }

	    DynamicMatrix(size_t M, size_t N, Type const& val)
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
    };

    template <typename Type, bool SO>
    struct EigenType<DynamicMatrix<Type, SO>>
    {
        typedef typename DynamicMatrix<Type, SO>::Base type;
    };

}