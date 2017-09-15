#pragma once

#define USE_BLAZE 0

#if USE_BLAZE
#include <blaze/Math.h>

namespace tmpc {


// Import blaze classed into tmpc namespace.
using namespace blaze;

}	// namespace tmpc

#else

#include "matrix/EigenAdapter.hpp"

#endif
