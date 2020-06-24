#pragma once


namespace tmpc
{
    template <
        typename Key,
        typename CustomVector,
        typename PtrPropertyMap,
        typename SizePropertyMap
    >
    class VectorPtrPropertyMap
    {
    public:
        VectorPtrPropertyMap(PtrPropertyMap ptr_map, SizePropertyMap size_map)
        :   ptrMap_{ptr_map}
        ,   sizeMap_{size_map}
        {
        }


        template <typename T>
        friend void put(VectorPtrPropertyMap const& pm, Key k, T const& val)
        {
            CustomVector v(get(pm.ptrMap_, k), get(pm.sizeMap_, k));
            v = val;
        }


        friend auto get(VectorPtrPropertyMap const& pm, Key k)
        {
            return CustomVector(get(pm.ptrMap_, k), get(pm.sizeMap_, k));
        }


    private:
        PtrPropertyMap ptrMap_;
        SizePropertyMap sizeMap_;
    };


    template <
        typename Key,
        typename CustomVector,
        typename PtrPropertyMap,
        typename SizePropertyMap
    >
    inline auto makeVectorPtrPropertyMap(PtrPropertyMap index_map, SizePropertyMap size_map)
    {
        return VectorPtrPropertyMap<Key, CustomVector, PtrPropertyMap, SizePropertyMap>(index_map, size_map);
    }
}