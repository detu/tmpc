#pragma once

#include <cmath>
#include <optional>


namespace tmpc
{
    /// @brief Unwrap data representing a cyclic quantity --
    /// for example, angle from an angle sensor.
    template <typename Real>
    class Unwrap
    {
    public:
        Unwrap(Real period = Real {2. * M_PI})
        :   period_(period)
        ,   halfPeriod_(period / 2.)
        {
        }


        template <typename T>
        Real operator()(T const& val)
        {
            if (last_)
            {
                Real const diff = val - *last_;
                Real const n = std::floor((diff + halfPeriod_) / period_);                
                *last_ += diff - period_ * n;
            }
            else
            {
                last_ = val;
            }

            return *last_;
        }


        Real period() const
        {
            return period_;
        }


        Real halfPeriod() const
        {
            return halfPeriod_;
        }


    private:
        Real const period_;
        Real const halfPeriod_;
        std::optional<Real> last_;
    };
}