#pragma once

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

    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, and object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		QpStageBase() = default;
		QpStageBase(QpStageBase const&) = default;
    };
}