#pragma once

#define USE_BLAZE 0

#if USE_BLAZE
#include <blaze/Math.h>

namespace tmpc 
{
    // Import blaze classed into tmpc namespace.
    using namespace blaze;
}

#else

#include "matrix/EigenAdaptor.hpp"

namespace tmpc 
{
    // Import blaze classed into tmpc namespace.
    using namespace eigen_adaptor;
}

#endif
