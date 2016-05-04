/*
 * qpDUNESSolution.hpp
 *
 *  Created on: May 4, 2016
 *      Author: kotlyar
 */

#ifndef INCLUDE_QP_QPDUNESSOLUTION_HPP_
#define INCLUDE_QP_QPDUNESSOLUTION_HPP_

#include <Eigen/Dense>

#include <vector>

namespace camels
{
	//
	// Manages data layout for qpDUNES solver output.
	//
	template<unsigned NX_, unsigned NU_>
	class qpDUNESSolution
	{
	public:
		typedef unsigned size_type;

		static size_type const NX = NX_;
		static size_type const NU = NU_;
		static size_type const NZ = NX + NU;

		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NZ, 1> StateInputVector;

		qpDUNESSolution(size_type nt)
		:	_data(NZ * nt + NX)
		,	_nt(nt)
		{
		}

		Eigen::Map<StateInputVector> w(unsigned i)
		{
			if (!(i < _nt))
				throw std::out_of_range("qpDUNESSolution<>::w(): index is out of range");

			return Eigen::Map<StateInputVector>(_data.data() + i * NZ);
		}

		Eigen::Map<StateInputVector const> w(unsigned i) const
		{
			if (!(i < _nt))
				throw std::out_of_range("qpDUNESSolution<>::w(): index is out of range");

			return Eigen::Map<StateInputVector const>(_data.data() + i * NZ);
		}

		Eigen::Map<StateVector> x(unsigned i)
		{
			if (!(i < _nt + 1))
				throw std::out_of_range("qpDUNESSolution<>::x(): index is out of range");

			return Eigen::Map<StateVector>(_data.data() + i * NZ);
		}

		Eigen::Map<StateVector const> x(unsigned i) const
		{
			if (!(i < _nt + 1))
				throw std::out_of_range("qpDUNESSolution<>::x(): index is out of range");

			return Eigen::Map<StateVector const>(_data.data() + i * NZ);
		}

		Eigen::Map<StateVector> wend()
		{
			return Eigen::Map<StateVector>(_data.data() + _nt * NZ);
		}

		Eigen::Map<StateVector const> wend() const
		{
			return Eigen::Map<StateVector const>(_data.data() + _nt * NZ);
		}

		void shift()
		{
			std::copy_n(_data.begin() + NZ, (_nt - 1) * NZ + NX, _data.begin());
		}

		qpDUNESSolution<NX_, NU_>& operator+=(qpDUNESSolution<NX_, NU_> const& rhs)
		{
			if (rhs.nT() != nT())
				throw std::invalid_argument("qpDUNESSolution<>::operator+=(): arguments have different sizes!");

			std::transform(_data.cbegin(), _data.cend(), rhs._data.cbegin(), _data.begin(), std::plus<double>());

			return *this;
		}

		size_type const nX() const noexcept { return NX; }
		size_type const nU() const noexcept { return NU; }
		size_type const nT() const noexcept { return _nt; }

		// qpDUNES interface
		//
		// Return pointer to data layout as expected by qpDUNES.
		double * data() { return _data.data(); }

	private:
		size_type const _nt;

		// _data stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _data;
	};
}

#endif /* INCLUDE_QP_QPDUNESSOLUTION_HPP_ */