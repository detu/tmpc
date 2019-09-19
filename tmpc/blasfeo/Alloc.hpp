#pragma once

#include <blasfeo_stdlib.h>

#include <tmpc/SizeT.hpp>


namespace tmpc :: blasfeo
{
    inline void * malloc_align(size_t size)
    {
        void * ptr = nullptr;
        blasfeo_malloc_align(&ptr, size);

        return ptr;
    }


    inline void free_align(void * ptr)
    {
        blasfeo_free_align(ptr);
    }
}