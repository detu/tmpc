#pragma once

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor 
{
    template <typename MT, typename OP>
    struct ForEachExpr 
    :   Eigen::MatrixBase<ForEachExpr<MT, OP>>
    {
        ForEachExpr(MT const& arg, OP const& op)
        :   arg_(arg)
        ,   op_(op)
        {
        }

        typedef typename Eigen::internal::ref_selector<ForEachExpr>::type Nested;

        typedef Eigen::Index Index;
        Index rows() const { return arg_.rows(); }
        Index cols() const { return arg_.rows(); }

        typedef typename Eigen::internal::ref_selector<MT>::type ArgTypeNested;

        ArgTypeNested arg_;
        OP const& op_;
    };
}

template <typename MT, typename OP>
struct ::Eigen::internal::traits<tmpc::eigen_adaptor::ForEachExpr<MT, OP>>
{
    typedef Eigen::Dense StorageKind;
    typedef Eigen::MatrixXpr XprKind;
    typedef typename MT::StorageIndex StorageIndex;
    typedef typename MT::Scalar Scalar;
    enum { 
        Flags = MT::Flags,
        RowsAtCompileTime = MT::RowsAtCompileTime,
        ColsAtCompileTime = MT::RowsAtCompileTime,
        MaxRowsAtCompileTime = MT::MaxRowsAtCompileTime,
        MaxColsAtCompileTime = MT::MaxRowsAtCompileTime,
        InnerStrideAtCompileTime = MT::InnerStrideAtCompileTime,
        OuterStrideAtCompileTime = MT::OuterStrideAtCompileTime
    };
};

template <typename MT, typename OP>
struct ::Eigen::internal::evaluator<tmpc::eigen_adaptor::ForEachExpr<MT, OP>>
:   evaluator_base<tmpc::eigen_adaptor::ForEachExpr<MT, OP>>
{
    typedef tmpc::eigen_adaptor::ForEachExpr<MT, OP> XprType;
    typedef typename nested_eval<MT, XprType::ColsAtCompileTime>::type ArgTypeNested;
    typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
    typedef typename XprType::CoeffReturnType CoeffReturnType;
    enum { 
        CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
        Flags = MT::Flags
    };
    
    evaluator(const XprType& xpr)
        : argImpl_(xpr.arg_)
    { }

    CoeffReturnType coeff(Index row, Index col) const
    {
        return argImpl_.op_(argImpl_.arg_(row, col));
    }

    evaluator<ArgTypeNestedCleaned> argImpl_;
};
