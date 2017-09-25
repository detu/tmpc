#pragma once

#include "ExpressionBase.hpp"
#include "Eval.hpp"

namespace tmpc
{
    template <typename OP, typename A>
    class UnaryOp
    :   public ExpressionBase<UnaryOp<OP, A>>
    {
    public:
        using Arg = A;
        using Op = OP;

        constexpr UnaryOp(Arg const& arg)
        :   arg_(arg)
        {        
        }

        auto const& arg() const
        {
            return arg_;
        }

    private:
        Arg arg_;
    };


    template <typename OP, typename A, class... Types>
    decltype(auto) eval(UnaryOp<OP, A> const& expr, std::tuple<Types...> const& arg)
    {
        return OP::eval(eval(expr.arg(), arg));
    }
}