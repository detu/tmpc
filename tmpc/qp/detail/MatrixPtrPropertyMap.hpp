#pragma once


namespace tmpc :: detail
{
    template <
        typename Key,
        typename CustomMatrix,
        typename PtrPropertyMap,
        typename SizePropertyMap
    >
    class MatrixPtrPropertyMap
    {
    public:
        MatrixPtrPropertyMap(PtrPropertyMap ptr_map, SizePropertyMap size_map)
        :   ptrMap_{ptr_map}
        ,   sizeMap_{size_map}
        {
        }


        template <typename T>
        friend void put(MatrixPtrPropertyMap const& pm, Key k, T const& val)
        {
            auto const sz = get(pm.sizeMap_, k);
            CustomMatrix m(get(pm.ptrMap_, k), sz.first, sz.second);
            m = val;
        }


        friend auto get(MatrixPtrPropertyMap const& pm, Key k)
        {
            auto const sz = get(pm.sizeMap_, k);
            return CustomMatrix(get(pm.ptrMap_, k), sz.first, sz.second);
        }


    private:
        PtrPropertyMap ptrMap_;
        SizePropertyMap sizeMap_;
    };


    template <
        typename Key,
        typename CustomMatrix,
        typename PtrPropertyMap,
        typename SizePropertyMap
    >
    inline auto makeMatrixPtrPropertyMap(PtrPropertyMap ptr_map, SizePropertyMap size_map)
    {
        return MatrixPtrPropertyMap<Key, CustomMatrix, PtrPropertyMap, SizePropertyMap>(ptr_map, size_map);
    }
}