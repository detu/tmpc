#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <iostream>
#include <cmath>

#include <tmpc/Matrix.hpp>

template <class Base>
class MatrixPrintWrap : public Base {
    friend void PrintTo(const MatrixPrintWrap &m, ::std::ostream *o) {
        *o << "\n" << m;
    }
};

/*
 * Makes the Eigen3 matrix classes printable from GTest checking macros like EXPECT_EQ.
 * Usage: EXPECT_EQ(print_wrap(a), print_wrap(b))
 * Taken from this post: http://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class Base>
const MatrixPrintWrap<Base> &print_wrap(const Base &base) {
    return static_cast<const MatrixPrintWrap<Base> &>(base);
}

MATCHER_P(FloatNearPointwise, tol, "Out of range") {
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
	MatrixApproxEquality(double tolerance)
	:	tolerance_(tolerance)
	{
	}

	template <typename MatrixA, typename MatrixB>
	bool operator()(MatrixA const& lhs, MatrixB const& rhs)
	{
		return l1Norm(lhs - rhs) <= tolerance_;
	}

private:
	double tolerance_;
};
