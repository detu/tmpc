#pragma once

#include <memory>


namespace tmpc :: blasfeo
{
    /// \brief Alignment for BLASFEO data arrays
    inline static size_t constexpr alignment()
    {
        return 0x40;
    }


    namespace detail
    {
        template <typename T>
        T * alignedAlloc(std::size_t bytes)
        {
            return reinterpret_cast<T *>(std::aligned_alloc(alignment(), bytes));
        }


        /// \brief Deleter for aligned data arrays
        struct AlignedDeleter
        {
            void operator()(void * data)
            {
                std::free(data);
            }
        };
    }
}