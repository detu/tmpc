#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc
{
	template <typename Real_>
	class TrajectoryPoint
	{
	public:
		typedef Real_ Real;

		decltype(auto) x()
		{
			return full(x_);
		}

		DynamicVector<Real> const& x() const
		{
			return x_;
		}

		decltype(auto) u()
		{
			return full(u_);
		}

		DynamicVector<Real> const& u() const
		{
			return u_;
		}

		decltype(auto) w()
		{
			return full(w_);
		}

		DynamicVector<Real> const& w() const
		{
			return w_;
		}

		decltype(auto) y()
		{
			return full(y_);
		}

		DynamicVector<Real> const& y() const
		{
			return y_;
		}

	private:
		DynamicVector<Real> x_;
		DynamicVector<Real> u_;
		DynamicVector<Real> w_;
		DynamicVector<Real> y_;
	};

	template <typename Real>
	inline std::ostream& operator<<(std::ostream& os, TrajectoryPoint<Real> const& p)
	{
		os << "x=" << trans(p.x()) 
			<< "\tu=" << trans(p.u()) 
			<< "\tw=" << trans(p.w()) 
			<< "\ty=" << trans(p.y());
	}
}
