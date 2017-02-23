#pragma once

#include "Types.hpp"

#include <initializer_list>

namespace tmpc
{
    using std::initializer_list;

    template <typename Type>
    inline size_t determineColumns(initializer_list<initializer_list<Type>> list) noexcept
    {
        size_t cols( 0UL );
        for( const auto& rowList : list )
            cols = std::max( cols, rowList.size() );
        return cols;
    }
}