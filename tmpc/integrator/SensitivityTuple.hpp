#pragma once

#include <tmpc/integrator/SensitivityTupleBase.hpp>

namespace tmpc
{
    template <typename Kernel>
    class SensitivityTuple
    :   public SensitivityTupleBase<SensitivityTuple<Kernel>>
    {
    public:
        SensitivityTuple(size_t ny, size_t nx, size_t nu)
        :   value_ {ny}
        ,   sensX_ {ny, nx}
        ,   sensU_ {ny, nu}
        {            
        }

        auto const& implValue() const
        {
            return value_;
        }

        template <typename T>
        void implValue(T const& val)
        {
            value_ = val;
        }

        auto const& implSensX() const
        {
            return sensX_;
        }

        template <typename T>
        void implSensX(T const& val)
        {
            sensX_ = val;
        }

        auto const& implSensU() const
        {
            return sensU_;
        }

        template <typename T>
        void implSensU(T const& val)
        {
            sensU_ = val;
        }

    private:
        DynamicVector<Kernel> value_;
		DynamicMatrix<Kernel> sensX_;
		DynamicVector<Kernel> sensU_;
    };
}