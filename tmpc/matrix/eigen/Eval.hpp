#pragma once

#include "Eigen.hpp"

namespace Eigen
{
    // Analog of blaze::eval():
    //
    // "The eval function forces the evaluation of the given dense matrix expression dm. 
    // The function returns an expression representing the operation."
    //
    // For Eigen matrices just returns its argument.
    //
    // Normally, you want to use evaluate().
    template <typename T>
    decltype(auto) eval(DenseBase<T> const& expr)
    {
        return expr;
    }


    // Evaluates the given matrix expression.
    template <typename T>
    decltype(auto) evaluate(DenseBase<T> const& expr)
    {
        return expr.eval();
    }
}