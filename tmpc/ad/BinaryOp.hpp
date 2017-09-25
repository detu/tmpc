#pragma once

#include "ExpressionBase.hpp"
#include "Eval.hpp"

namespace tmpc
{
    template <typename OP, typename L, typename R>
    class BinaryOp
    :   public ExpressionBase<BinaryOp<OP, L, R>>
    {
    public:
        using Left = L;
        using Right = R;
        using Op = OP;

        constexpr BinaryOp(Left const& left, Right const& right)
        :   left_(left)
        ,   right_(right)
        {        
        }

        auto const& left() const
        {
            return left_;
        }

        auto const& right() const
        {
            return right_;
        }

        size_t implOperationCount() const
        {
            return tmpc::operationCount(left_) + tmpc::operationCount(right_) + 1;
        }
        
    private:
        Left left_;
        Right right_;
    };


    template <typename OP, typename L, typename R, class... Types>
    decltype(auto) eval(BinaryOp<OP, L, R> const& expr, std::tuple<Types...> const& arg)
    {
        return OP::eval(eval(expr.left(), arg), eval(expr.right(), arg));
    }


    template <typename OP, typename L, typename R, std::size_t N>
    auto constexpr diff(BinaryOp<OP, L, R> const& expr, Variable<N> const& u)
    {
        // Chain rule: d_OP(left, right)/d_left * d_left/d_u + d_OP(left, right)/d_right * d_right/d_u
        return OP::diffL(expr.left(), expr.right()) * diff(expr.left(), u)
            + OP::diffR(expr.left(), expr.right()) * diff(expr.right(), u);
    }
}