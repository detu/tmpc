#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/qp/OcpQpExpressionBase.hpp>
#include <tmpc/Matrix.hpp>

#include <boost/range/iterator_range_core.hpp>

#include <stdexcept>


namespace tmpc
{
    ///
	/// \brief CRTP base for OCP QP vertex properties.
	///
	template <typename Derived>
	class OcpQpVertexBase
	{
	public:
		// -----------------------------------------------------------
		// Cost
		// -----------------------------------------------------------
        decltype(auto) Q() const 
		{ 
			return derived().impl_Q(); 
		}

		
		template <typename T> 
		Derived& Q(const T& q) 
		{ 
			derived().impl_Q(q);
			return derived();
		}


		template <typename T>
		Derived& Q(std::initializer_list<std::initializer_list<T>> val) 
		{
			derived().impl_Q(val);
			return derived();
		}

		
		decltype(auto) R() const 
		{
			return derived().impl_R(); 
		}


		template <typename T> 
		Derived& R(const T& r) 
		{ 
			derived().impl_R(r);
			return derived();
		}


		template <typename T>
		Derived& R(std::initializer_list<std::initializer_list<T>> val) 
		{
			derived().impl_R(val);
			return derived();
		}


		decltype(auto) S() const 
		{	
			return derived().impl_S(); 
		}


		template <typename T> 
		Derived& S(const T& s) 
		{ 
			derived().impl_S(s); 
			return derived(); 
		}


		template <typename T>
		Derived& S(std::initializer_list<std::initializer_list<T>> val) 
		{
			derived().impl_S(val);
			return derived();
		}


		decltype(auto) q() const 
		{	
			return derived().impl_q(); 
		}


		template <typename T> 
		Derived& q(const T& q) 
		{ 
			derived().impl_q(q); 
			return derived(); 
		}


		template <typename T>
		Derived& q(std::initializer_list<T> val) 
		{
			derived().impl_q(val); 
			return derived(); 
		}


		decltype(auto) r() const 
		{	
			return derived().impl_r(); 
		}


		template <typename T> 
		Derived& r(const T& r) 
		{ 
			derived().impl_r(r); 
			return derived(); 
		}


		template <typename T>
		Derived& r(std::initializer_list<T> val) 
		{
			derived().impl_r(val); 
			return derived(); 
		}


		// -----------------------------------------------------------
		// Soft constraints cost
		// -----------------------------------------------------------
		decltype(auto) Zl() const { return derived().impl_Zl(); }
		template <typename T> Derived& Zl(const T& val) { derived().impl_Zl(val); return derived(); }

		decltype(auto) Zu() const {	return derived().impl_Zu(); }
		template <typename T> Derived& Zu(const T& val) { derived().impl_Zu(val); return derived(); }

		decltype(auto) zl() const { return derived().impl_zl(); }
		template <typename T> Derived& zl(const T& val) { derived().impl_zl(val); return derived(); }

		decltype(auto) zu() const {	return derived().impl_zu(); }
		template <typename T> Derived& zu(const T& val) { derived().impl_zu(val); return derived(); }

		// -----------------------------------------------------------
		// Linear constraints
		// -----------------------------------------------------------
		decltype(auto) C() const {	return derived().C(); }
		template <typename T> Derived& C(const T& c) { derived().C(c); return derived(); }

		decltype(auto) D() const {	return derived().D(); }
		template <typename T> Derived& D(const T& d) { derived().D(d); return derived(); }

		decltype(auto) lbd() const { return derived().lbd(); }
		template <typename T> Derived& lbd(const T& lbd) { derived().lbd(lbd); return derived(); }

		decltype(auto) ubd() const { return derived().ubd(); }
		template <typename T> Derived& ubd(const T& ubd) { derived().ubd(ubd); return derived(); }

		// -----------------------------------------------------------
		// Soft linear constraints
		// -----------------------------------------------------------
		decltype(auto) idxs() const 
		{ 
			return derived().impl_idxs(); 
		}
		
		
		template <typename T> 
		Derived& idxs(T const& val) 
		{
			if (val.size() != derived().impl_idxs().size())
				throw std::invalid_argument("Soft constraints index size does not match");

			derived().impl_idxs(val);
			return derived();
		}


		template <typename T> 
		Derived& idxs(std::initializer_list<T> val) 
		{ 
			derived().impl_idxs(val);
			return derived();
		}

		// -----------------------------------------------------------
		// Bound constraints
		// -----------------------------------------------------------
		decltype(auto) lbu() const { return derived().lbu(); }
		template <typename T> Derived& lbu(const T& lbu) { derived().lbu(lbu); return derived(); }

		decltype(auto) ubu() const { return derived().ubu(); }
		template <typename T> Derived& ubu(const T& ubu) { derived().ubu(ubu); return derived(); }

		decltype(auto) lbx() const { return derived().lbx(); }
		template <typename T> Derived& lbx(const T& lbx) { derived().lbx(lbx); return derived(); }		

		decltype(auto) ubx() const { return derived().ubx(); }
        template <typename T> Derived& ubx(const T& ubx) { derived().ubx(ubx); return derived(); }


		/// Problem size
        decltype(auto) size() const
        {
            return derived().size();
		}
		
		
		template <typename Other>
		OcpQpVertexBase& operator=(OcpQpVertexBase<Other> const& rhs)
		{
			Q(rhs.Q());
			R(rhs.R());
			S(rhs.S());
			q(rhs.q());
			r(rhs.r());
			Zl(rhs.Zl());
			Zu(rhs.Zu());
			zl(rhs.zl());
			zu(rhs.zu());
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


		/*
		/// \brief Expression assignment.
		template <typename Expr>
		OcpQpVertexBase& operator=(OcpQpExpressionBase<Expr> const& expr)
		{
			expr.evalTo(*this);
		}
		*/


		// Set all data to sNaN
		Derived& setNaN()
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
			C(sNaN<Real>());
			D(sNaN<Real>());
			lbd(sNaN<Real>());
			ubd(sNaN<Real>());
			lbu(sNaN<Real>());
			ubu(sNaN<Real>());
			lbx(sNaN<Real>());
			ubx(sNaN<Real>());

			return derived();
		}
		
		// Set the Gauss-Newton approximation of the hessian Hessian and the gradient.
		template <typename ResidualVector, typename CMatrix, typename DMatrix>
		Derived& gaussNewtonCostApproximation(ResidualVector const& res, CMatrix const& C, DMatrix const& D)
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

			return derived();
		}

		
		// Set upper and lower input bounds relative to a point.
		template <typename Vector1, typename Vector2, typename Vector3>
		Derived& relativeInputBounds(Vector1 const& u, Vector2 const& lu, Vector3 const& uu)
		{
			lbu(lu - u);
			ubu(uu - u);

			return derived();
		}

		// Set upper and lower state bounds relative to a point.
		template <typename Vector1, typename Vector2, typename Vector3>
		Derived& relativeStateBounds(Vector1 const& x, Vector2 const& lx, Vector3 const& ux)
		{
			lbx(lx - x);
			ubx(ux - x);

			return derived();
		}


		// Set path constraint inequality relative to a point.
		template <typename Vector1, typename Vector2, typename Matrix1, typename Matrix2, typename Vector3, typename Vector4>
		Derived& relativePathConstraints(Vector1 const& x, Vector2 const& u, 
			Matrix1 const& C, Matrix2 const& D, Vector3 const& ld, Vector4 const& ud)
		{
			this->C(C);
			this->D(D);
			this->lbd(ld - C * x - D * u);
			this->ubd(ud - C * x - D * u);

			return derived();
		}


		// Set upper and lower input bounds.
		template <typename Vector1, typename Vector2>
		Derived& inputBounds(Vector1 const& lu, Vector2 const& uu)
		{
			lbu(lu);
			ubu(uu);

			return derived();
		}

		// Set upper and lower state bounds.
		template <typename Vector1, typename Vector2>
		Derived& stateBounds(Vector1 const& lx, Vector2 const& ux)
		{
			lbx(lx);
			ubx(ux);

			return derived();
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename IteratorRange, typename Matrix, typename Vector>
		Derived& softConstraints(IteratorRange const& idxs, Matrix const& Z, Vector const& z)
		{
			this->idxs(idxs);
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

			return derived();
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename Matrix, typename Vector>
		Derived& slackPenalty(Matrix const& Z, Vector const& z)
		{
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

			return derived();
		}


		/// \brief Set soft constraints.
		///
		/// Sets slack index and sets slack Hessian and gradient as Zl = Zu = Z, zl = zu = z.
		template <typename T, typename Matrix, typename Vector>
		Derived& softConstraints(std::initializer_list<T> idxs, Matrix const& Z, Vector const& z)
		{
			this->idxs(idxs);
			Zl(Z);
			Zu(Z);
			zl(z);
			zu(z);

			return derived();
		}

    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		OcpQpVertexBase() = default;
		OcpQpVertexBase(OcpQpVertexBase const&) = default;

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
    inline void randomize(OcpQpVertexBase<QP>& qp)
    {
        using Kernel = typename QP::Kernel;
		using DynamicMatrix = DynamicMatrix<Kernel>;
		using DynamicVector = DynamicVector<Kernel>;
		typename Kernel::template Rand<DynamicMatrix> rand_matrix;
		typename Kernel::template Rand<DynamicVector> rand_vector;

		auto const sz = qp.size();

		{
			DynamicMatrix H = rand_matrix.generate(sz.nx() + sz.nu(), sz.nx() + sz.nu());
			H *= ctrans(H);

			qp.Q(submatrix(H, 0, 0, sz.nx(), sz.nx()));
			qp.R(submatrix(H, sz.nx(), sz.nx(), sz.nu(), sz.nu()));
			qp.S(submatrix(H, 0, sz.nx(), sz.nx(), sz.nu()));
		}

		qp.q(rand_vector.generate(sz.nx()));
		qp.r(rand_vector.generate(sz.nu()));
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