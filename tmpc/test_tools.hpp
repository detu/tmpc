#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <iostream>
#include <cmath>

#include <tmpc/Matrix.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
	namespace detail
	{
		template <typename T>
		class ForcePrint 
		: 	public T 
		{
			friend void PrintTo(const ForcePrint &m, ::std::ostream *o) 
			{
				*o << "\n" << m;
			}
		};
	}


	/*
	* Makes the Eigen3 matrix classes printable from GTest checking macros like EXPECT_EQ.
	* Usage: EXPECT_EQ(forcePrint(a), forcePrint(b))
	* Taken from this post: http://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
	*/
	template <typename T>
	decltype(auto) forcePrint(T const &base) 
	{
		return static_cast<detail::ForcePrint<T> const&>(base);
	}


	MATCHER_P(FloatNearPointwise, tol, "Out of range") 
	{
		return (std::get<0>(arg) > std::get<1>(arg) - tol && std::get<0>(arg) < std::get<1>(arg) + tol) ;
	}


	// Provides standard iterators over a matrix
	template <typename Matrix>
	class MatrixContainerAdaptor
	{
	public:
		typedef size_t size_type;
		typedef typename Matrix::ElementType value_type;
		typedef value_type const& const_reference;

		MatrixContainerAdaptor(Matrix const& m) : m_(m) {};

		class const_iterator
		{
		public:
			const_iterator(Matrix const &m, size_type i, size_type j) : i_(i), j_(j), m_(m) {}

			const_iterator& operator++()
			{
				if (++i_ >= m_.rows())
				{
					i_ = 0;
					++j_;
				}

				return *this;
			}

			decltype(auto) operator*() const
			{
				return m_(i_, j_);
			}

		private:
			size_type i_;
			size_type j_;
			Matrix const& m_;
		};

		const_iterator begin() const { return const_iterator(m_, 0,         0); }
		const_iterator end  () const { return const_iterator(m_, 0, m_.cols()); }
		size_type size() const { return m_.rows() * m_.cols(); }

		Matrix const& getMatrix() const { return m_; }

	private:
		Matrix const& m_;
	};


	template <typename Matrix>
	std::ostream& operator<<(std::ostream& os, MatrixContainerAdaptor<Matrix> const& ma)
	{
		return os << ma.getMatrix();
	}


	template <class Matrix>
	MatrixContainerAdaptor<Matrix> as_container(Matrix const& m)
	{
		return MatrixContainerAdaptor<Matrix>(m);
	}


	class MatrixApproxEquality
	{
	public:
		MatrixApproxEquality(double abs_tol, double rel_tol = 0.)
		:	absTol_(abs_tol)
		,	relTol_(rel_tol)
		{
		}


		template <typename MT1, bool SO1, typename MT2, bool SO2>
		bool operator()(blaze::Matrix<MT1, SO1> const& lhs, blaze::Matrix<MT2, SO2> const& rhs) const
		{
			size_t const M = rows(lhs);
			size_t const N = columns(lhs);

			if (rows(rhs) != M || columns(rhs) != N)
				throw std::invalid_argument("Matrix size mismatch in MatrixApproxEquality");

			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					if (abs((~lhs)(i, j) - (~rhs)(i, j)) > absTol_ + relTol_ * abs((~rhs)(i, j)))
						return false;

			return true;
		}


		template <typename VT1, typename VT2, bool TF>
		bool operator()(blaze::Vector<VT1, TF> const& lhs, blaze::Vector<VT2, TF> const& rhs) const
		{
			size_t const N = size(lhs);

			if (size(rhs) != N)
				throw std::invalid_argument("Vector size mismatch in MatrixApproxEquality");

			for (size_t j = 0; j < N; ++j)
				if (abs((~lhs)[j] - (~rhs)[j]) > absTol_ + relTol_ * abs((~rhs)[j]))
					return false;

			return true;
		}


	private:
		double const absTol_;
		double const relTol_;
	};
}