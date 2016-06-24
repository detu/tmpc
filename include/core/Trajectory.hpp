#pragma once

#include <Eigen/Dense>

#include <stdexcept>
#include <vector>

namespace tmpc
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

	private:
		Derived& derived() { return static_cast<Derived&>(*this); }
		Derived const& derived() const { return static_cast<Derived const&>(*this); }
	};

	template <unsigned NX, unsigned NU>
	class Trajectory
	{
	public:
		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NU, 1> InputVector;

		Trajectory(std::size_t n) : _x(n + 1), _u(n) {}

		template <typename StateVector_, typename InputVector_>
		Trajectory(std::size_t n, Eigen::MatrixBase<StateVector_> const& x, Eigen::MatrixBase<InputVector_> const& u)
		: _x(n + 1, x), _u(n, u) {}

		StateVector const& get_x(std::size_t i) const { return _x.at(i); }
		template <typename Matrix> void set_x(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { _x.at(i) = val; }

		InputVector const& get_u(std::size_t i) const { return _u.at(i); }
		template <typename Matrix> void set_u(std::size_t i, Eigen::MatrixBase<Matrix> const& val) { _u.at(i) = val; }

		std::size_t nT() const { return _u.size(); }

		template<typename Other>
		Trajectory& operator+=(Other const& rhs)
		{
			if (rhs.nT() != nT())
				throw std::invalid_argument("Trajectory::operator+=(): arguments must have same number of time steps!");

			for (std::size_t i = 0; i < nT() + 1; ++i)
				_x[i] += rhs.get_x(i);

			for (std::size_t i = 0; i < nT(); ++i)
			    _u[i] += rhs.get_u(i);

			return *this;
		}

	private:
		std::vector<StateVector> _x;
		std::vector<InputVector> _u;
	};

	template<unsigned NX_, unsigned NU_>
	void shift(Trajectory<NX_, NU_>& tr)
	{
		for (std::size_t i = 0; i < tr.nT(); ++i)
			tr.set_x(i, tr.get_x(i + 1));

		for (std::size_t i = 0; i + 1 < tr.nT(); ++i)
			tr.set_u(i, tr.get_u(i + 1));
	}

	// The following gives me "unable to deduce template arguments" error:
	template<unsigned NX_, unsigned NU_>
	inline std::ostream& operator<<(std::ostream& os, Trajectory<NX_, NU_> const& tr)
	{
		for (std::size_t i = 0; i < tr.nT(); ++i)
			os << tr.get_x(i).transpose() << "\t" << tr.get_u(i).transpose() << std::endl;

		return os << tr.get_x(tr.nT()).transpose() << std::endl;
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
