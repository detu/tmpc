/*** 
 * This is a temporary solution to determine the Kernel type and the Real type
 * used by, for example, a QP solver. We need to get rid of it once we got rid of kernels.
*/

#pragma once


namespace tmpc
{
    /// @brief Get Kernel type of Whatever
    template <typename Whatever>
    struct KernelOf
    {
        using type = typename Whatever::Kernel;
    };


    /// @brief Get Real type of Whatever
    template <typename Whatever>
    struct RealOf
    {
        using type = typename Whatever::Real;
    };
}