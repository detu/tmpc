#pragma once

#include <tmpc/Matrix.hpp>

#include <stdexcept>
#include <vector>

namespace tmpc
{
	template <unsigned NX, unsigned NU, unsigned NW = 0, unsigned NY = 0>
	class Trajectory
	{
	public:
		typedef StaticVector<double, NX> StateVector;
		typedef StaticVector<double, NU> InputVector;
		typedef StaticVector<double, NW> DisturbanceVector;
		typedef StaticVector<double, NY> MeasurementVector;

		/**
		 * \brief Default constructor.
		 */
		Trajectory(std::size_t n) : x_(n + 1), u_(n), w_(n), y_(n) {}

		/**
		 * \brief Constructor for MPC.
		 */
		template <typename StateVector_, typename InputVector_>
		Trajectory(std::size_t n, StateVector_ const& x, InputVector_ const& u)
		: x_(n + 1, x), u_(n, u), w_(n), y_(n) {}

		/**
		 * \brief Constructor for MHE.
		 */
		template <typename StateVector_, typename InputVector_, typename DisturbanceVector_, typename MeasurementVector_>
		Trajectory(std::size_t n, StateVector_ const& x, InputVector_ const& u,	DisturbanceVector_ const& w, MeasurementVector_ const& y)
		: x_(n + 1, x), u_(n, u), w_(n, w), y_(n, y) {}

		StateVector const& get_x(std::size_t i) const { return x_.at(i); }
		template <typename Matrix> void set_x(std::size_t i, Matrix const& val) { x_.at(i) = val; }

		/**
		 * \brief Get input value at stage k
		 */
		InputVector const& get_u(std::size_t k) const
		{
			return u_.at(k);
		}

		/**
		 * \brief Set input value at stage k
		 */
		template <typename Vector>
		void set_u(std::size_t k, Vector const& val)
		{
			u_.at(k) = val;
		}

		/**
		 * \brief Set part of input value at stage k starting from element i
		 */
		template <typename Vector>
		void set_u(std::size_t k, unsigned i, Vector const& val)
		{
			subvector(u_.at(k), i, val.size()) = val;
		}

		DisturbanceVector const& get_w(std::size_t i) const { return w_.at(i); }
		template <typename Matrix> void set_w(std::size_t i, Matrix const& val) { w_.at(i) = val; }

		MeasurementVector const& get_y(std::size_t i) const { return y_.at(i); }
		template <typename Matrix> void set_y(std::size_t i, Matrix const& val) { y_.at(i) = val; }

		std::size_t nT() const { return u_.size(); }

	private:
		std::vector<StateVector      > x_;
		std::vector<InputVector      > u_;
		std::vector<DisturbanceVector> w_;
		std::vector<MeasurementVector> y_;
	};

	template<unsigned NX_, unsigned NU_, unsigned NW_, unsigned NY_>
	void shift(Trajectory<NX_, NU_, NW_, NY_>& tr)
	{
		for (std::size_t i = 0; i < tr.nT(); ++i)
			tr.set_x(i, tr.get_x(i + 1));

		for (std::size_t i = 0; i + 1 < tr.nT(); ++i)
		{
			tr.set_u(i, tr.get_u(i + 1));
			tr.set_w(i, tr.get_w(i + 1));
			tr.set_y(i, tr.get_y(i + 1));
		}
	}

	// The following gives me "unable to deduce template arguments" error:
	template<unsigned NX_, unsigned NU_, unsigned NW_, unsigned NY_>
	inline std::ostream& operator<<(std::ostream& os, Trajectory<NX_, NU_, NW_, NY_> const& tr)
	{
		for (std::size_t i = 0; i < tr.nT(); ++i)
			os << tr.get_x(i).transpose() << "\t" << tr.get_u(i).transpose()
			   << tr.get_w(i).transpose() << "\t" << tr.get_y(i).transpose() << std::endl;

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
