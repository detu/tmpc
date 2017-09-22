#pragma once

#include <tmpc/Math.hpp>

namespace tmpc
{
    ///
	/// \brief CRTP base for QP stage classes.
	///
	template <typename Derived>
	class QpStageBase
	{
	public:
		Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}

		Derived const& derived() const
		{
			return static_cast<Derived const&>(*this);	
        }
        
        decltype(auto) Q() const { return derived().Q(); }
		template <typename T> void Q(const T& q) { derived().Q(q); }

		decltype(auto) R() const {	return derived().R(); }
		template <typename T> void R(const T& r) { derived().R(r); }

		decltype(auto) S() const {	return derived().S(); }
		template <typename T> void S(const T& s) { derived().S(s); }

		decltype(auto) q() const {	return derived().q(); }
		template <typename T> void q(const T& q) { derived().q(q); }

		decltype(auto) r() const {	return derived().r(); }
		template <typename T> void r(const T& r) { derived().r(r); }

		decltype(auto) A() const {	return derived().A(); }
		template <typename T> void A(const T& a) { derived().A(a); }

		decltype(auto) B() const {	return derived().B(); }
		template <typename T> void B(const T& b) { derived().B(b); }

		decltype(auto) b() const { return derived().b(); }
		template <typename T> void b(const T& b) { derived().b(b); }

		decltype(auto) C() const {	return derived().C(); }
		template <typename T> void C(const T& c) { derived().C(c); }

		decltype(auto) D() const {	return derived().D(); }
		template <typename T> void D(const T& d) { derived().D(d); }

		decltype(auto) lbd() const { return derived().lbd(); }
		template <typename T> void lbd(const T& lbd) { derived().lbd(lbd); }

		decltype(auto) ubd() const { return derived().ubd(); }
		template <typename T> void ubd(const T& ubd) { derived().ubd(ubd); }

		decltype(auto) lbu() const { return derived().lbu(); }
		template <typename T> void lbu(const T& lbu) { derived().lbu(lbu); }

		decltype(auto) ubu() const { return derived().ubu(); }
		template <typename T> void ubu(const T& ubu) { derived().ubu(ubu); }

		decltype(auto) lbx() const { return derived().lbx(); }
		template <typename T> void lbx(const T& lbx) { derived().lbx(lbx); }		

		decltype(auto) ubx() const { return derived().ubx(); }
        template <typename T> void ubx(const T& ubx) { derived().ubx(ubx); }
        
        decltype(auto) size() const
        {
            return derived().size();
		}

		// Set all data to sNaN
		void setNaN()
		{
			using Kernel = typename Derived::Kernel;
			using Real = typename Kernel::Real;
	
			Q(sNaN<Real>());
			R(sNaN<Real>());
			S(sNaN<Real>());
			q(sNaN<Real>());
			r(sNaN<Real>());
			A(sNaN<Real>());
			B(sNaN<Real>());
			b(sNaN<Real>());
			C(sNaN<Real>());
			D(sNaN<Real>());
			lbd(sNaN<Real>());
			ubd(sNaN<Real>());
			lbu(sNaN<Real>());
			ubu(sNaN<Real>());
			lbx(sNaN<Real>());
			ubx(sNaN<Real>());
		}
		
		// Set the Gauss-Newton approximation of the hessian Hessian and the gradient.
		template <typename ResidualVector, typename CMatrix, typename DMatrix>
		void gaussNewtonCostApproximation(ResidualVector const& res, CMatrix const& C, DMatrix const& D)
		{
			// H = G^T G
			//   = [Q S
			//      S R]
			//

			Q(trans(C) * C);
			R(trans(D) * D);
			S(trans(C) * D);

			// g = 2 * (y_bar - y_hat)^T * W * G
			// g = [q; r]
			q(trans(C) * res);
			r(trans(D) * res);
		}

		// Set A, B and b to represent a linearized shooting equality
		// of the form \Delta x_{k+1} = \frac{\dif f}{\dif x}(x_k,u_k)\, \Delta x_{k} 
		//	+ \frac{\dif f}{\dif u}(x_k,u_k)\, \Delta u_{k} + f(x_{k},u_{k}) - x_{k+1}
		template <typename Vector1, typename Matrix1, typename Matrix2, typename Vector2>
		void linearizedShootingEquality(Vector1 const& f, Matrix1 const& Jx, Matrix2 const& Ju, Vector2 const& x_plus)
		{
			A(Jx);
			B(Ju);
			b(f - x_plus);
		}

		// Set upper and lower bounds relative to a point.
		template <typename Vector1, typename Vector2, typename Vector3, typename Vector4, typename Vector5, typename Vector6>
		void relativeBounds(Vector1 const& x, Vector2 const& u, Vector3 const& lx, Vector4 const& lu, Vector5 const& ux, Vector6 const& uu)
		{
			lbx(lx - x);
			ubx(ux - x);
			lbu(lu - u);
			ubu(uu - u);
		}

		// Set upper and lower bounds.
		template <typename Vector1, typename Vector2, typename Vector3, typename Vector4>
		void bounds(Vector1 const& lx, Vector2 const& lu, Vector3 const& ux, Vector4 const& uu)
		{
			lbx(lx);
			ubx(ux);
			lbu(lu);
			ubu(uu);
		}

    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, and object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		QpStageBase() = default;
		QpStageBase(QpStageBase const&) = default;
    };
}