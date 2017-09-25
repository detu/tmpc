#pragma once

#include "ExpressionBase.hpp"
//#include "Eval.hpp"
//#include "Diff.hpp"
#include "Zero.hpp"

/*
namespace tmpc
{
    template <typename T>
    class Value
    :   public ExpressionBase<Value<T>>
    {
    public:
        Value(T const& val)
        :   value_ {val}
        {            
        }

        auto const& value() const
        {
            return value_;
        }

    private:
        T value_;
    };
    
    
    template <typename T, class... Types>
    decltype(auto) eval(Value<T> const& expr, std::tuple<Types...> const& arg)
    {
        return expr.value();
    }


    template <typename T, std::size_t N2>
    decltype(auto) constexpr diff(Value<T> const& v, Variable<N2> const& u)
    {
        return Zero {};
    }
}
*/