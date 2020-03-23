#pragma once

#include "QpSolverException.hpp"

#include <algorithm>
#include <exception>


namespace tmpc 
{
	/**
	 * \brief CRTP base for QP workspace classes.
	 */
	template <typename Derived>
	class QpWorkspaceBase
	{
	public:
		void solve()
		{
			derived().impl_solve();
		}


		decltype(auto) solution() const
		{
			return derived().impl_solution();	
		}


		std::string solverName() const
		{
			return derived().impl_solverName();
		}
		
	protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		QpWorkspaceBase() = default;
		QpWorkspaceBase(QpWorkspaceBase const&) = default;

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
}
