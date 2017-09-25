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

    
    template <typename OP, typename A, std::size_t N, typename S, class... Types>
    decltype(auto) constexpr diff(UnaryOp<OP, A> const& expr, Variable<N> const& u, std::tuple<Types...> const& arg, S const& sens)
    {
        return OP::diff(eval(expr.arg(), arg), sens);
    }
}