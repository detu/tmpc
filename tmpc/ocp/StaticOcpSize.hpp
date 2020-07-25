#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/Exception.hpp>


namespace tmpc 
{
	/**
	 * \brief Defines sizes of a statically-sized OCP.
	*/
	template <size_t NX, size_t NU, size_t NC = 0>
	class StaticOcpSize
	{
	public:
		explicit StaticOcpSize(OcpTree const& graph)
		:	graph_ {graph}
		{
		}


		/**
		 * \brief The OCP graph
		 */
		auto const& graph() const noexcept
		{
			return graph_;
		}


		/**
		* \brief Number of states.
		*/
		size_t nx(OcpVertex v) const
		{
			return NX;
		}


		/**
		* \brief Number of inputs.
		*/
		size_t nu(OcpVertex v) const
		{
			return out_degree(v, graph_) > 0 ? NU : 0;
		}


		/**
		* \brief Number of constraints.
		*/
		size_t nc(OcpVertex v) const
		{
			return NC;
		}


		/**
		* \brief Number of soft constraints.
		*/
		size_t ns(OcpVertex v) const
		{
			return 0;
		}


	private:
		OcpTree const graph_;
	};

	
	template <size_t NX1, size_t NU1, size_t NC1, size_t NX2, size_t NU2, size_t NC2>
	inline bool operator==(StaticOcpSize<NX1, NU1, NC1> const& a, StaticOcpSize<NX2, NU2, NC2> const& b)
	{
		return NX1 == NX2 && NU1 == NU2 && NC1 == NC2 && a.graph() == b.graph();
	}

	
	template <size_t NX1, size_t NU1, size_t NC1, size_t NX2, size_t NU2, size_t NC2>
	inline bool operator!=(StaticOcpSize<NX1, NU1, NC1> const& a, StaticOcpSize<NX2, NU2, NC2> const& b)
	{
		return !(a == b);
	}


	/// @brief Total number of variables in an OCP
	///
	template <size_t NX, size_t NU, size_t NC>
	inline size_t numVariables(StaticOcpSize<NX, NU, NC> const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nx(v) + s.nu(v) + 2 * s.ns(v);

		return n;
	}


	/// @brief Total number of state variables in an OCP.
	///
	template <size_t NX, size_t NU, size_t NC>
	inline size_t numStates(StaticOcpSize<NX, NU, NC> const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nx(v);

		return n;
	}


	/// @brief Total number of control variables in an OCP.
	///
	template <size_t NX, size_t NU, size_t NC>
	inline size_t numControls(StaticOcpSize<NX, NU, NC> const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nu(v);

		return n;
	}


	/// @brief Total number of equality constraints in an OCP.
	///
	template <size_t NX, size_t NU, size_t NC>
	inline size_t numEqualities(StaticOcpSize<NX, NU, NC> const& s)
	{
		size_t n = 0;
		for (auto e : edges(s.graph()))
			n += s.nx(target(e, s.graph()));

		return n;
	}


	/**
	* \brief Total number of rows in all inequalities like lbg[k] <= C[k]*x[k] + D[k]*u[k] <= ubg[k]
	*/
	template <size_t NX, size_t NU, size_t NC>
	inline size_t numInequalities(StaticOcpSize<NX, NU, NC> const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nc(v);

		return n;
	}	
	

	/**
	 * Resulting OCP size after full condensing.
	 * 
	 * TODO: move to Condensing.hpp
	 */
	template <size_t NX, size_t NU, size_t NC>
	inline DynamicOcpSize condensedSize(StaticOcpSize<NX, NU, NC> const& s)
	{
		auto const nx0 = s.nx(root(s.graph()));

		return DynamicOcpSize {{
			nx0,
			numControls(s),
			numInequalities(s) + numStates(s) - nx0
			// TODO: handle the number of soft constraints properly.
		}};
	}
}
