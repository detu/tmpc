#pragma once

#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>


namespace tmpc
{
    using boost::iterator_property_map;
    using boost::function_property_map;
    using boost::associative_property_map;
    using boost::const_associative_property_map;
    using boost::transform_value_property_map;
    
    using boost::make_function_property_map;
    using boost::make_iterator_property_map;
    using boost::make_transform_value_property_map;
    
    using boost::read_write_property_map_tag;
    using boost::readable_property_map_tag;
    using boost::read_write_property_map_tag;

    using boost::property_traits;
    
    
    /// @brief Copy values from property map src to property map dst for the keys defined by the iterator range keys.
    template <typename PropMapSrc, typename PropMapDst, typename KeyRange>
    inline void copyProperty(PropMapSrc src, PropMapDst dst, KeyRange keys)
    {
        for (auto key : keys)
            put(dst, key, get(src, key));
    }
}
