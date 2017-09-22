#pragma once


namespace tmpc
{
	///
	/// \brief CRTP base for QP stage solution classes.
	///
	template <typename Derived>
	class OcpSolutionBase
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

		decltype(auto) x() const { return derived().x();	}
		decltype(auto) u() const { return derived().u();	}
		decltype(auto) pi() const	{ return derived().pi(); }
		decltype(auto) lam_lbu() const { return derived().lam_lbu(); }
		decltype(auto) lam_ubu() const { return derived().lam_ubu(); }
		decltype(auto) lam_lbx() const { return derived().lam_lbx(); }
		decltype(auto) lam_ubx() const { return derived().lam_ubx(); }
		decltype(auto) lam_lbd() const { return derived().lam_lbd(); }
		decltype(auto) lam_ubd() const { return derived().lam_ubd(); }

	protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, and object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		OcpSolutionBase() = default;
		OcpSolutionBase(OcpSolutionBase const&) = default;
	};
}
