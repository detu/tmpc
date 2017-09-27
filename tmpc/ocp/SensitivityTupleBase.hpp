#pragma once

namespace tmpc
{
    template <typename Derived>
    class SensitivityTupleBase
    {
    public:
        decltype(auto) value() const
        {
            return derived().implValue();
        }

        template <typename T>
        void value(T const& val) const
        {
            return derived().implValue(val);
        }

        decltype(auto) sensX() const
        {
            return derived().implSensX();
        }

        template <typename T>
        void sensX(T const& val) const
        {
            return derived().implSensX(val);
        }

        decltype(auto) sensU() const
        {
            return derived().implSensU();
        }

        template <typename T>
        void sensU(T const& val) const
        {
            return derived().implSensU(val);
        }

    private:
        Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}

		Derived const& derived() const
		{
			return static_cast<Derived const&>(*this);	
		}
    };
}