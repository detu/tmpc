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
            pm.impl_put(k, val);
        }


        friend decltype(auto) get(VectorPropertyMap const& pm, Key k)
        {
            return pm.impl_get(k);
        }


    private:
        void impl_put(Key k, Vector const& val) const
        {
            auto const sz = get(sizeMap_, k);
            if (sz != val.size())
                throw std::invalid_argument("Invalid vector size");

            setter_(val.data(), get(indexMap_, k));
        }


        auto impl_get(Key k) const
        {
            auto const sz = get(sizeMap_, k);

            Vector val(sz);
            getter_(val.data(), get(indexMap_, k));

            return val;
        }


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