#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/qp/OcpQpExpressionBase.hpp>
#include <tmpc/Matrix.hpp>

#include <boost/range/iterator_range_core.hpp>

#include <stdexcept>


namespace tmpc
{
    ///
	/// \brief CRTP base for OCP QP stage classes.
	///
	template <typename Derived>
	class OcpQpEdgeBase
	{
	public:
		// -----------------------------------------------------------
		// Shooting equalities
		// -----------------------------------------------------------
		decltype(auto) A() const {	return derived().A(); }
		template <typename T> Derived& A(const T& a) { derived().A(a); return derived(); }

		decltype(auto) B() const {	return derived().B(); }
		template <typename T> Derived& B(const T& b) { derived().B(b); return derived(); }

		decltype(auto) b() const { return derived().b(); }
		template <typename T> Derived& b(const T& b) { derived().b(b); return derived(); }

		/// Problem size
        decltype(auto) size() const
        {
            return derived().size();
		}

		template <typename Other>
		OcpQpEdgeBase& operator=(OcpQpEdgeBase<Other> const& rhs)
		{
			A(rhs.A());
			B(rhs.B());
			b(rhs.b());
			
			return *this;
		}


		/*
		/// \brief Expression assignment.
		template <typename Expr>
		OcpQpEdgeBase& operator=(OcpQpExpressionBase<Expr> const& expr)
		{
			expr.evalTo(*this);
		}
		*/


		// Set all data to sNaN
		Derived& setNaN()
		{
			using Kernel = typename Derived::Kernel;
			using Real = typename Kernel::Real;

			A(sNaN<Real>());
			B(sNaN<Real>());
			b(sNaN<Real>());

			return derived();
		}
		
		
		// Set A, B and b to represent a linearized shooting equality
		// of the form \Delta x_{k+1} = \frac{\dif f}{\dif x}(x_k,u_k)\, \Delta x_{k} 
		//	+ \frac{\dif f}{\dif u}(x_k,u_k)\, \Delta u_{k} + f(x_{k},u_{k}) - x_{k+1}
		template <typename Vector1, typename Matrix1, typename Matrix2, typename Vector2>
		Derived& linearizedShootingEquality(Vector1 const& f, Matrix1 const& Jx, Matrix2 const& Ju, Vector2 const& x_plus)
		{
			A(Jx);
			B(Ju);
			b(f - x_plus);

			return derived();
		}


		// Set A, B and b to specified values.
		template <typename Matrix1, typename Matrix2, typename Vector>
		Derived& shootingEquality(Matrix1 const& A, Matrix2 const& B, Vector const& b)
		{
			this->A(A);
			this->B(B);
			this->b(b);

			return derived();
		}


    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		OcpQpEdgeBase() = default;
		OcpQpEdgeBase(OcpQpEdgeBase const&) = default;

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
    inline void randomize(OcpQpEdgeBase<QP>& qp)
    {
        using Kernel = typename QP::Kernel;
		using DynamicMatrix = DynamicMatrix<Kernel>;
		using DynamicVector = DynamicVector<Kernel>;
		typename Kernel::template Rand<DynamicMatrix> rand_matrix;
		typename Kernel::template Rand<DynamicVector> rand_vector;

		auto const sz = qp.size();
		auto const nx_next = rows(qp.A());

		qp.A(rand_matrix.generate(nx_next, sz.nx()));
		qp.B(rand_matrix.generate(nx_next, sz.nu()));
		qp.b(rand_vector.generate(nx_next));
    }
}