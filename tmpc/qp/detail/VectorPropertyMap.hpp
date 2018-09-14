#pragma once

#include <stdexcept>


namespace tmpc :: detail
{
    template <
        typename Key,
        typename Vector,
        typename IndexPropertyMap,
        typename SizePropertyMap,
        typename SetterFunc,
        typename GetterFunc
    >
    class VectorPropertyMap
    {
    public:
        VectorPropertyMap(IndexPropertyMap index_map, SizePropertyMap size_map, SetterFunc setter, GetterFunc getter)
        :   indexMap_{index_map}
        ,   sizeMap_{size_map}
        ,   setter_{setter}
        ,   getter_{getter}
        {
        }


        friend void put(VectorPropertyMap const& pm, Key k, Vector const& val)
        {
            auto const sz = get(pm.sizeMap_, k);
            if (sz != val.size())
                throw std::invalid_argument("Invalid vector size");

            pm.setter_(val.data(), get(pm.indexMap_, k));
        }


        friend decltype(auto) get(VectorPropertyMap const& pm, Key k)
        {
            auto const sz = get(pm.sizeMap_, k);

            Vector val(sz);
            pm.getter_(val.data(), get(pm.indexMap_, k));

            return val;
        }


    private:
        IndexPropertyMap indexMap_;
        SizePropertyMap sizeMap_;
        SetterFunc setter_;
        GetterFunc getter_;
    };


    template <
        typename Key,
        typename Vector,
        typename IndexPropertyMap,
        typename SizePropertyMap,
        typename SetterFunc,
        typename GetterFunc
    >
    inline auto makeVectorPropertyMap(IndexPropertyMap index_map, SizePropertyMap size_map, SetterFunc setter, GetterFunc getter)
    {
        return VectorPropertyMap<Key, Vector, IndexPropertyMap, SizePropertyMap, SetterFunc, GetterFunc>(index_map, size_map, setter, getter);
    }
}