#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/qp/OcpQpExpressionBase.hpp>
#include <tmpc/qp/OcpQpVertexBase.hpp>
#include <tmpc/qp/OcpQpEdgeBase.hpp>
#include <tmpc/Matrix.hpp>

#include <boost/range/iterator_range_core.hpp>

#include <stdexcept>


namespace tmpc
{
    ///
	/// \brief CRTP base for OCP QP stage classes.
	///
	template <typename Derived>
	class OcpQpBase
	:	public OcpQpVertexBase<Derived>
	,	public OcpQpEdgeBase<Derived>
	{
	public:
		/// Problem size
        decltype(auto) size() const
        {
            return derived().size();
		}


		template <typename Other>
		OcpQpBase& operator=(OcpQpBase<Other> const& rhs)
		{
			if (size() != rhs.size())
				throw std::invalid_argument("OcpQp assignment size mismatch");

			static_cast<OcpQpVertexBase<Derived>&>(*this) = rhs;
			static_cast<OcpQpEdgeBase<Derived>&>(*this) = rhs;
			
			return *this;
		}


		/// \brief Expression assignment.
		template <typename Expr>
		OcpQpBase& operator=(OcpQpExpressionBase<Expr> const& expr)
		{
			expr.evalTo(*this);
		}


		// Set all data to sNaN
		OcpQpBase& setNaN()
		{
			using Kernel = typename Derived::Kernel;
			using Real = typename Kernel::Real;

			static_cast<OcpQpVertexBase<Derived>&>(*this).setNaN();
			static_cast<OcpQpEdgeBase<Derived>&>(*this).setNaN();

			return *this;
		}


    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		OcpQpBase() = default;
		OcpQpBase(OcpQpBase const&) = default;


	private:
		Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}


		Derived const& derived() const
		{
			return static_cast<Derived const&>(*this);	
        }
    };


	///
    /// \brief Randomize a QP.
    ///
    template <typename QP>
    inline void randomize(OcpQpBase<QP>& qp)
    {
        using Kernel = typename QP::Kernel;
		using DynamicMatrix = DynamicMatrix<Kernel>;
		using DynamicVector = DynamicVector<Kernel>;
		typename Kernel::template Rand<DynamicMatrix> rand_matrix;
		typename Kernel::template Rand<DynamicVector> rand_vector;

		randomize(static_cast<OcpQpVertexBase<QP>&>(qp));
		randomize(static_cast<OcpQpEdgeBase<QP>&>(qp));
    }
}