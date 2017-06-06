#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc
{
	template <typename Scalar_>
	class TrajectoryPoint
	{
	public:
		typedef Scalar_ Scalar;

		decltype(auto) x()
		{
			return full(x_);
		}

		DynamicVector<Scalar> const& x() const
		{
			return x_;
		}

		decltype(auto) u()
		{
			return full(u_);
		}

		DynamicVector<Scalar> const& u() const
		{
			return u_;
		}

		decltype(auto) w()
		{
			return full(w_);
		}

		DynamicVector<Scalar> const& w() const
		{
			return w_;
		}

		decltype(auto) y()
		{
			return full(y_);
		}

		DynamicVector<Scalar> const& y() const
		{
			return y_;
		}

	private:
		DynamicVector<Scalar> x_;
		DynamicVector<Scalar> u_;
		DynamicVector<Scalar> w_;
		DynamicVector<Scalar> y_;
	};

	template <typename Scalar>
	inline std::ostream& operator<<(std::ostream& os, TrajectoryPoint<Scalar> const& p)
	{
		os << "x=" << trans(p.x()) 
			<< "\tu=" << trans(p.u()) 
			<< "\tw=" << trans(p.w()) 
			<< "\ty=" << trans(p.y());
	}
}
