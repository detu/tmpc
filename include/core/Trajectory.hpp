#pragma once

namespace camels
{
	template<typename Derived_, unsigned NX_, unsigned NU_>
	class TrajectoryBase
	{
	public:
		typedef unsigned size_type;

		static size_type const NX = NX_;
		static size_type const NU = NU_;
		static size_type const NZ = NX + NU;
	};

	template<typename Trajectory, typename StateVector, typename InputVector>
	Trajectory ConstantTrajectory(typename Trajectory::size_type nt, StateVector const& x, InputVector const& u)
	{
		Trajectory tr(nt);

		for (typename Trajectory::size_type i = 0; i < nt; ++i)
		{
			tr.x(i) = x;
			tr.u(i) = u;
		}

		tr.wend() = x;

		return tr;
	}
}
