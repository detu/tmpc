#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc
{
	template <typename Scalar_>
	class MpcTrajectoryPoint
	{
	public:
		typedef Scalar_ Scalar;

		MpcTrajectoryPoint(size_t NX, size_t NU)
		:	x_(NX)
		,	u_(NU)
		{
		}

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

	private:
		DynamicVector<Scalar> x_;
		DynamicVector<Scalar> u_;
	};

	template <typename Scalar>
	inline std::ostream& operator<<(std::ostream& os, MpcTrajectoryPoint<Scalar> const& p)
	{
		os << "x=" << trans(p.x()) << "\tu=" << trans(p.u());
	}
}
