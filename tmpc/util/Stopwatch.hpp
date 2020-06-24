#pragma once

#include <chrono>


namespace tmpc
{
    /// @brief Utility class to measure time intervals.
    ///
    template <typename Clock_ = std::chrono::steady_clock>
    class Stopwatch
    {
    public:
        /// @brief Clock type used by the Stopwatch
        using Clock = Clock_;

        /// @brief Duration type used by the Stopwatch, synonym to Clock::duration.
        using Duration = typename Clock::duration;


        /// @brief Constructor
        ///
        /// Records the current time as start time.
        /// Remembers the reference to the dur variable.
        Stopwatch(Duration& dur)
        :   duration_ {dur}
        ,   startTime_ {Clock::now()}
        {
            static_assert(Clock::is_steady);
        }


        /// @brief Destructor
        ///
        /// Subtracts the previously recorded start time from the current time, 
        /// and stores the result in the Duration variable who's reference was previously passed to the constructor.
        ~Stopwatch()
        {
            duration_ = Clock::now() - startTime_;
        }
        

    private:
        Duration& duration_;
        typename Clock::time_point startTime_;
    };
}