#pragma once

#include <stdexcept>

namespace camels
{
	template<typename Derived_, unsigned NX_, unsigned NU_>
	class TrajectoryBase
	{
	public:
		typedef unsigned size_type;
		typedef Derived_ Derived;

		/*
		typedef typename Derived_::StateVector StateVector;
		typedef typename Derived_::StateInputVector StateInputVector;
		*/

		/*
		static size_type const NX = NX_;
		static size_type const NU = NU_;
		static size_type const NZ = NX + NU;
		*/

		decltype(auto) w(size_type i) const noexcept { return derived().w(i); }
		decltype(auto) w(size_type i) noexcept { return derived().w(i); }
		decltype(auto) wend() const noexcept { return derived().wend(); }
		decltype(auto) wend() noexcept { return derived().wend(); }
		decltype(auto) x(size_type i) const noexcept { return derived().x(i); }
		decltype(auto) x(size_type i) noexcept { return derived().x(i); }

		size_type nT() const noexcept { return derived().nT(); }

		template<typename OtherDerived_>
		Derived& operator+=(TrajectoryBase<OtherDerived_, NX_, NU_> const& rhs)
		{
			if (rhs.nT() != nT())
				throw std::invalid_argument("TrajectoryBase<>::operator+=(): arguments have different sizes!");

			for (size_type i = 0; i + 1 < nT(); ++i)
				w(i) += rhs.w(i);

			wend() += rhs.wend();

			return derived();
		}

	private:
		Derived& derived() { return static_cast<Derived&>(*this); }
		Derived const& derived() const { return static_cast<Derived const&>(*this); }
	};

	// Creates trajectory with nt steps. At every step, control is set equal to u and state is set equal to x.
	//
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

	template<typename Trajectory, unsigned NX_, unsigned NU_>
	void shift(TrajectoryBase<Trajectory, NX_, NU_>& tr)
	{
		for (typename Trajectory::size_type i = 0; i + 1 < tr.nT(); ++i)
			tr.w(i) = tr.w(i + 1);

		tr.x(tr.nT() - 1) = tr.wend();
	}

	/*
	template<typename Trajectory1_, typename Trajectory2_, unsigned NX_, unsigned NU_>
	Trajectory1_& operator+=(qpDUNESSolution<NX_, NU_> const& rhs)
	{
		if (rhs.nT() != nT())
			throw std::invalid_argument("qpDUNESSolution<>::operator+=(): arguments have different sizes!");

		std::transform(_data.cbegin(), _data.cend(), rhs._data.cbegin(), _data.begin(), std::plus<double>());

		return *this;
	}
	*/
}
