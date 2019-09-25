#pragma once

#include <blasfeo_stdlib.h>
#include <blasfeo_d_aux_ext_dep.h>

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


    inline void allocate_mat(size_t m, size_t n, blasfeo_dmat& mat)
    {
        blasfeo_allocate_dmat(m, n, &mat);
    }


    inline void free_mat(blasfeo_dmat& mat)
    {
        blasfeo_free_dmat(&mat);
    }
}