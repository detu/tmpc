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
	class OcpQpBase
	{
	public:
		// -----------------------------------------------------------
		// Cost
		// -----------------------------------------------------------
        decltype(auto) Q() const { return derived().Q(); }
		
		template <typename T> 
		OcpQpBase& Q(const T& q) 
		{ 
			derived().Q(q);
			return *this;
		}

		decltype(auto) R() const {	return derived().R(); }
		template <typename T> OcpQpBase& R(const T& r) 
		{ 
			derived().R(r);
			return *this;
		}

		decltype(auto) S() const {	return derived().S(); }
		template <typename T> OcpQpBase& S(const T& s) { derived().S(s); return *this; }

		decltype(auto) q() const {	return derived().q(); }
		template <typename T> OcpQpBase& q(const T& q) { derived().q(q); return *this; }

		decltype(auto) r() const {	return derived().r(); }
		template <typename T> OcpQpBase& r(const T& r) { derived().r(r); return *this; }

		// -----------------------------------------------------------
		// Soft constraints cost
		// -----------------------------------------------------------
		decltype(auto) Zl() const { return derived().impl_Zl(); }
		template <typename T> OcpQpBase& Zl(const T& val) { derived().impl_Zl(val); return *this; }

		decltype(auto) Zu() const {	return derived().impl_Zu(); }
		template <typename T> OcpQpBase& Zu(const T& val) { derived().impl_Zu(val); return *this; }

		decltype(auto) zl() const { return derived().impl_zl(); }
		template <typename T> OcpQpBase& zl(const T& val) { derived().impl_zl(val); return *this; }

		decltype(auto) zu() const {	return derived().impl_zu(); }
		template <typename T> OcpQpBase& zu(const T& val) { derived().impl_zu(val); return *this; }

		// -----------------------------------------------------------
		// Shooting equalities
		// -----------------------------------------------------------
		decltype(auto) A() const {	return derived().A(); }
		template <typename T> OcpQpBase& A(const T& a) { derived().A(a); return *this; }

		decltype(auto) B() const {	return derived().B(); }
		template <typename T> OcpQpBase& B(const T& b) { derived().B(b); return *this; }

		decltype(auto) b() const { return derived().b(); }
		template <typename T> OcpQpBase& b(const T& b) { derived().b(b); return *this; }

		// -----------------------------------------------------------
		// Linear constraints
		// -----------------------------------------------------------
		decltype(auto) C() const {	return derived().C(); }
		template <typename T> OcpQpBase& C(const T& c) { derived().C(c); return *this; }

		decltype(auto) D() const {	return derived().D(); }
		template <typename T> OcpQpBase& D(const T& d) { derived().D(d); return *this; }

		decltype(auto) lbd() const { return derived().lbd(); }
		template <typename T> OcpQpBase& lbd(const T& lbd) { derived().lbd(lbd); return *this; }

		decltype(auto) ubd() const { return derived().ubd(); }
		template <typename T> OcpQpBase& ubd(const T& ubd) { derived().ubd(ubd); return *this; }

		// -----------------------------------------------------------
		// Soft linear constraints
		// -----------------------------------------------------------
		decltype(auto) idxs() const 
		{ 
			return derived().impl_idxs(); 
		}
		
		
		template <typename T> 
		OcpQpBase& idxs(T const& val) 
		{
			if (val.size() != size().ns())
				throw std::invalid_argument("Soft constraints index size does not match");

			derived().impl_idxs(val);
			return *this;
		}


		template <typename T> 
		OcpQpBase& idxs(std::initializer_list<T> val) 
		{ 
			derived().impl_idxs(val);
			return *this;
		}

		// -----------------------------------------------------------
		// Bound constraints
		// -----------------------------------------------------------
		decltype(auto) lbu() const { return derived().lbu(); }
		template <typename T> OcpQpBase& lbu(const T& lbu) { derived().lbu(lbu); return *this; }

		decltype(auto) ubu() const { return derived().ubu(); }
		template <typename T> OcpQpBase& ubu(const T& ubu) { derived().ubu(ubu); return *this; }

		decltype(auto) lbx() const { return derived().lbx(); }
		template <typename T> OcpQpBase& lbx(const T& lbx) { derived().lbx(lbx); return *this; }		

		decltype(auto) ubx() const { return derived().ubx(); }
        template <typename T> OcpQpBase& ubx(const T& ubx) { derived().ubx(ubx); return *this; }
		
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

			Q(rhs.Q());
			R(rhs.R());
			S(rhs.S());
			q(rhs.q());
			r(rhs.r());
			Zl(rhs.Zl());
			Zu(rhs.Zu());
			zl(rhs.zl());
			zu(rhs.zu());
			A(rhs.A());
			B(rhs.B());
			b(rhs.b());
			C(rhs.C());
			D(rhs.D());
			lbx(rhs.lbx());
			ubx(rhs.ubx());
			lbu(rhs.lbu());
			ubu(rhs.ubu());
			lbd(rhs.lbd());
			ubd(rhs.ubd());
			idxs(rhs.idxs());

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
	
			Q(sNaN<Real>());
			R(sNaN<Real>());
			S(sNaN<Real>());
			q(sNaN<Real>());
			r(sNaN<Real>());
			Zl(sNaN<Real>());
			Zu(sNaN<Real>());
			zl(sNaN<Real>());
			zu(sNaN<Real>());
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

			return *this;
		}
		
		// Set the Gauss-Newton approximation of the hessian Hessian and the gradient.
		template <typename ResidualVector, typename CMatrix, typename DMatrix>
		OcpQpBase& gaussNewtonCostApproximation(ResidualVector const& res, CMatrix const& C, DMatrix const& D)
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

			return *this;
		}

		// Set A, B and b to represent a linearized shooting equality
		// of the form \Delta x_{k+1} = \frac{\dif f}{\dif x}(x_k,u_k)\, \Delta x_{k} 
		//	+ \frac{\dif f}{\dif u}(x_k,u_k)\, \Delta u_{k} + f(x_{k},u_{k}) - x_{k+1}
		template <typename Vector1, typename Matrix1, typename Matrix2, typename Vector2>
		OcpQpBase& linearizedShootingEquality(Vector1 const& f, Matrix1 const& Jx, Matrix2 const& Ju, Vector2 const& x_plus)
		{
			A(Jx);
			B(Ju);
			b(f - x_plus);

			return *this;
		}


		// Set A, B and b to specified values.
		template <typename Matrix1, typename Matrix2, typename Vector>
		OcpQpBase& shootingEquality(Matrix1 const& A, Matrix2 const& B, Vector const& b)
		{
			this->A(A);
			this->B(B);
			this->b(b);

			return *this;
		}


		// Set upper and lower input bounds relative to a point.
		template <typename Vector1, typename Vector2, typename Vector3>
		OcpQpBase& relativeInputBounds(Vector1 const& u, Vector2 const& lu, Vector3 const& uu)
		{
			lbu(lu - u);
			ubu(uu - u);

			return *this;
		}

		// Set upper and lower state bounds relative to a point.
		template <typename Vector1, typename Vector2, typename Vector3>
		OcpQpBase& relativeStateBounds(Vector1 const& x, Vector2 const& lx, Vector3 const& ux)
		{
			lbx(lx - x);
			ubx(ux - x);

			return *this;
		}

		// Set upper and lower input bounds.
		template <typename Vector1, typename Vector2>
		OcpQpBase& inputBounds(Vector1 const& lu, Vector2 const& uu)
		{
			lbu(lu);
			ubu(uu);

			return *this;
		}

		// Set upper and lower state bounds.
		template <typename Vector1, typename Vector2>
		OcpQpBase& stateBounds(Vector1 const& lx, Vector2 const& ux)
		{
			lbx(lx);
			ubx(ux);

			return *this;
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename IteratorRange, typename Matrix, typename Vector>
		OcpQpBase& softConstraints(IteratorRange const& idxs, Matrix const& Z, Vector const& z)
		{
			this->idxs(idxs);
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

			return *this;
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename Matrix, typename Vector>
		OcpQpBase& slackPenalty(Matrix const& Z, Vector const& z)
		{
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

			return *this;
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack index and sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename T, typename Matrix, typename Vector>
		OcpQpBase& softConstraints(std::initializer_list<T> idxs, Matrix const& Z, Vector const& z)
		{
			this->idxs(idxs);
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

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

		auto const sz = qp.size();
		auto const nx_next = rows(qp.A());

		{
			DynamicMatrix H = rand_matrix.generate(sz.nx() + sz.nu(), sz.nx() + sz.nu());
			H *= ctrans(H);

			qp.Q(submatrix(H, 0, 0, sz.nx(), sz.nx()));
			qp.R(submatrix(H, sz.nx(), sz.nx(), sz.nu(), sz.nu()));
			qp.S(submatrix(H, 0, sz.nx(), sz.nx(), sz.nu()));
		}

		qp.q(rand_vector.generate(sz.nx()));
		qp.r(rand_vector.generate(sz.nu()));
		qp.A(rand_matrix.generate(nx_next, sz.nx()));
		qp.B(rand_matrix.generate(nx_next, sz.nu()));
		qp.b(rand_vector.generate(nx_next));
		qp.C(rand_matrix.generate(sz.ns(), sz.nx()));
		qp.D(rand_matrix.generate(sz.ns(), sz.nu()));

		{
			DynamicVector const lbd = rand_vector.generate(sz.ns());
			DynamicVector const ubd = rand_vector.generate(sz.ns());
			qp.lbd(min(lbd, ubd));
			qp.ubd(max(lbd, ubd));
		}

		{
			DynamicVector const lbx = rand_vector.generate(sz.nx());
			DynamicVector const ubx = rand_vector.generate(sz.nx());
			qp.lbx(min(lbx, ubx));
			qp.ubx(max(lbx, ubx));
		}

		{
			DynamicVector const lbu = rand_vector.generate(sz.nu());
			DynamicVector const ubu = rand_vector.generate(sz.nu());
			qp.lbu(min(lbu, ubu));
			qp.ubu(max(lbu, ubu));
		}

		{
			DynamicMatrix const Z = rand_matrix.generate(sz.ns(), sz.ns());
			qp.Zl(ctrans(Z) * Z);
		}

		{
			DynamicMatrix const Z = rand_matrix.generate(sz.ns(), sz.ns());
			qp.Zu(ctrans(Z) * Z);
		}
			
		qp.zl(rand_vector.generate(sz.ns()));
		qp.zu(rand_vector.generate(sz.ns()));
    }
}