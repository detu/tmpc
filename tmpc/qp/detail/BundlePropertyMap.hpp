#pragma once

#include <tmpc/core/PropertyMap.hpp>


namespace tmpc :: detail
{
    template <
        typename Value,
        typename BundleMap
    >
    class BundlePropertyMap
    {
    public:
        using key_type = typename property_traits<BundleMap>::key_type;
        using value_type = Value;
        using reference = Value&;
        using category = read_write_property_map_tag;
        using Bundle = typename property_traits<BundleMap>::value_type;


        BundlePropertyMap(Value Bundle:: * field, BundleMap bundle_map)
        :   field_{field}
        ,   bundleMap_{bundle_map}
        {
        }


        template <typename V>
        friend void put(BundlePropertyMap const& pm, key_type k, V const& val)
        {
            Value& ref = pm.bundleMap_[k].*pm.field_;
            
            if (shape(val) != shape(ref))
                throw std::invalid_argument("Invalid size of a vector or a matrix");

            ref = val;
        }


        friend auto& get(BundlePropertyMap const& pm, key_type k)
        {
            return pm.bundleMap_[k].*pm.field_;
        }


    private:
        Value Bundle:: * field_;
        BundleMap bundleMap_;
    };
}