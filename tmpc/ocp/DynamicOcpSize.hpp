#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/Exception.hpp>

#include <initializer_list>
#include <memory>
#include <algorithm>


namespace tmpc 
{
	struct OcpVertexSize
	{
		OcpVertexSize(size_t nx = 0, size_t nu = 0, size_t nc = 0, size_t ns = 0)
		:	nx_ {nx}
		,	nu_ {nu}
		,	nc_ {nc}
		,	ns_ {ns}
		{
		}


		size_t nx_ = 0;
		size_t nu_ = 0;
		size_t nc_ = 0;
		size_t ns_ = 0;
	};


	/**
	 * \brief Defines sizes of an OCP stage.
	 * 
	 * TODO: rename to DynamicOcpSize
	*/
	class DynamicOcpSize
	{
	public:
		/// @brief Construct size for a nominal OCP.
		//
		explicit DynamicOcpSize(std::initializer_list<OcpVertexSize> vertex_size);

		/// @brief Construct size for a tree-structured OCP.
		//
		explicit DynamicOcpSize(OcpTree const& graph, std::initializer_list<OcpVertexSize> vertex_size);


		/// @brief Construct size for a tree-structured OCP
		/// spcecifying sizes for each vertex individually.
		//
		template <typename SizeFunctor>
		explicit DynamicOcpSize(OcpTree const& graph, SizeFunctor sf)
		:	graph_ {graph}
		,	size_ {new OcpVertexSize[num_vertices(graph_)]}
		{
			for (auto v : vertices(graph_))
				size_[v] = sf(v, graph_);
		}


		/// @brief Construct OcpSize for a given tree 
		/// using constant nx, nu, nc for internal nodes.
		///
		explicit DynamicOcpSize(OcpTree const& graph, size_t nx, size_t nu,
			size_t nc = 0, size_t ns = 0, size_t nct = 0, bool first_state_empty = true);


		/// @brief Construct OcpSize for a single-scenario tree 
		/// using constant nx, nu, nc for internal nodes.
		///
		explicit DynamicOcpSize(size_t num_stages, size_t nx, size_t nu,
			size_t nc = 0, size_t ns = 0, size_t nct = 0, bool first_state_empty = true);


		/// @brief Copy from another \a OcpSize object
		///
		template <OcpSize Size>
		explicit DynamicOcpSize(Size const& size)
		:	graph_ {size.graph()}
		,	size_ {new OcpVertexSize[num_vertices(graph_)]}
		{
			for (auto v : vertices(graph_))
				size_[v] = OcpVertexSize {size.nx(v), size.nu(v), size.nc(v)};
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
			return size_[v].nx_;
		}


		/**
		* \brief Number of inputs.
		*/
		size_t nu(OcpVertex v) const
		{
			return size_[v].nu_;
		}


		/**
		* \brief Number of constraints.
		*/
		size_t nc(OcpVertex v) const
		{
			return size_[v].nc_;
		}


		/**
		* \brief Number of soft constraints.
		*/
		size_t ns(OcpVertex v) const
		{
			return size_[v].ns_;
		}


	private:
		OcpTree graph_;
		std::shared_ptr<OcpVertexSize[]> size_;
	};

	
	/// @brief Total number of variables in an OCP
	///
	inline size_t numVariables(DynamicOcpSize const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nx(v) + s.nu(v) + 2 * s.ns(v);

		return n;
	}


	/// @brief Total number of state variables in an OCP.
	///
	inline size_t numStates(DynamicOcpSize const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nx(v);

		return n;
	}


	/// @brief Total number of control variables in an OCP.
	///
	inline size_t numControls(DynamicOcpSize const& s)
	{
		size_t n = 0;
		for (auto v : vertices(s.graph()))
			n += s.nu(v);

		return n;
	}


	/// @brief Total number of equality constraints in an OCP.
	///
	inline size_t numEqualities(DynamicOcpSize const& s)
	{
		size_t n = 0;
		for (auto e : edges(s.graph()))
			n += s.nx(target(e, s.graph()));

		return n;
	}


	/**
	* \brief Total number of rows in all inequalities like lbg[k] <= C[k]*x[k] + D[k]*u[k] <= ubg[k]
	*/
	inline size_t numInequalities(DynamicOcpSize const& s)
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
	inline DynamicOcpSize condensedSize(DynamicOcpSize const& s)
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
