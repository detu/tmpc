#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc
{
	/**
	 * \brief A minimization variable of an OCP problem stage.
	 */
	template <typename Kernel>
	class OcpPoint
	{
	public:
		using Real = typename Kernel::Real;

		OcpPoint(size_t nx, size_t nu)
		:	x_ {nx}
		,	u_ {nu}
		{
		}

		template <typename VectorX, typename VectorU>
		OcpPoint(VectorX const& x, VectorU const& u)
		:	x_ {x}
		,	u_ {u}
		{
		}

		template <typename T>
		void x(T const& val)
		{
			noresize(x_) = val;
		}

		auto const& x() const
		{
			return x_;
		}

		template <typename T>
		void u(T const& val)
		{
			noresize(u_) = val;
		}

		auto const& u() const
		{
			return u_;
		}

	private:
		DynamicVector<Kernel> x_;
		DynamicVector<Kernel> u_;
	};

	template <typename Kernel>
	inline std::ostream& operator<<(std::ostream& os, OcpPoint<Kernel> const& p)
	{
		os << "x=" << trans(p.x()) << "\tu=" << trans(p.u());
	}
}
