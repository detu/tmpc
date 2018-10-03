#pragma once

#include <tmpc/Matrix.hpp>


namespace tmpc :: detail
{
    template <
        typename Key,
        typename Matrix,
        typename IndexPropertyMap,
        typename SizePropertyMap,
        typename SetterFunc,
        typename GetterFunc
    >
    class MatrixPropertyMap
    {
    public:
        MatrixPropertyMap(IndexPropertyMap index_map, SizePropertyMap size_map, SetterFunc setter, GetterFunc getter)
        :   indexMap_{index_map}
        ,   sizeMap_{size_map}
        ,   setter_{setter}
        ,   getter_{getter}
        {
        }


        friend void put(MatrixPropertyMap const& pm, Key k, Matrix const& val)
        {
            auto const sz = get(pm.sizeMap_, k);
            if (sz.first != rows(val) || sz.second != columns(val))
                throw std::invalid_argument("Invalid matrix size");

            pm.setter_(val.data(), spacing(val), get(pm.indexMap_, k));
        }


        friend decltype(auto) get(MatrixPropertyMap const& pm, Key k)
        {
            auto const sz = get(pm.sizeMap_, k);                
            Matrix m(sz.first, sz.second);

            pm.getter_(m.data(), spacing(m), get(pm.indexMap_, k));

            return m;
        }


    private:
        IndexPropertyMap indexMap_;
        SizePropertyMap sizeMap_;
        SetterFunc setter_;
        GetterFunc getter_;
    };


    template <
        typename Key,
        typename Matrix,
        typename IndexPropertyMap,
        typename SizePropertyMap,
        typename SetterFunc,
        typename GetterFunc
    >
    inline auto makeMatrixPropertyMap(IndexPropertyMap index_map, SizePropertyMap size_map, SetterFunc setter, GetterFunc getter)
    {
        return MatrixPropertyMap<Key, Matrix, IndexPropertyMap, SizePropertyMap, SetterFunc, GetterFunc>(index_map, size_map, setter, getter);
    }
}