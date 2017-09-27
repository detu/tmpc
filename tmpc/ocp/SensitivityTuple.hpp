#pragma once

#include <tmpc/ocp/SensitivityTupleBase.hpp>
#include <tmpc/Matrix.hpp>

#include <tuple>
#include <cassert>


namespace tmpc
{
    template <typename Kernel, size_t NV>
    class SensitivityTuple
    //:   public SensitivityTupleBase<SensitivityTuple<Kernel>>
    {
    public:
        /*
        template <typename ... Args>
        SensitivityTuple(size_t ny, Args ... nx)
        :   value_ (ny)
        ,   sens_ {DynamicMatrix<Kernel> {ny, nx} ...}
        {
        }
        */


        SensitivityTuple(size_t ny, std::initializer_list<size_t> const nx)
        :   value_ (ny)
        {
            // TODO: assert(nx.size() == NV); does not compile here; but why?
            assert(nx.size() == NV);

            auto nxx = nx.begin();
            for (auto& s : sens_)
                s.resize(ny, *nxx++);
        }


        auto const& value() const
        {
            return value_;
        }


        template <size_t N>
        auto const& sens() const
        {
            static_assert(N < NV);
            return sens_[N];
        }

        /*
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
        */

    private:
        DynamicVector<Kernel> value_;
		std::array<DynamicMatrix<Kernel>, NV> sens_;
    };
}