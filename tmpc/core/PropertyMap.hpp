#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>


namespace tmpc
{
    using boost::iterator_property_map;
    using boost::function_property_map;
    using boost::associative_property_map;
    
    using boost::make_function_property_map;
    using boost::transform_value_property_map;
    using boost::read_write_property_map_tag;
}
