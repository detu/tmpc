#pragma once

#include <tmpc/Math.hpp>

#include <tmpc/qp/OcpQp.hpp>

#include <vector>


namespace tmpc
{
	/// @brief Generates a random feasible OCP QP
	///
    template <OcpQp Qp>
	inline void randomize(Qp& qp)
	{
		using Real = typename Qp::Real;
		// Using matrix of the same storage order as in qp to workaround this issue:
		// https://bitbucket.org/blaze-lib/blaze/issues/364
		using DynamicMatrix = blaze::DynamicMatrix<Real, blaze::columnMajor>;
		using DynamicVector = blaze::DynamicVector<Real>;
		blaze::Rand<DynamicMatrix> rand_matrix;
		blaze::Rand<DynamicVector> rand_vector;
		
		auto const& g = qp.graph();
		auto const vert = vertices(g);

		// Record feasible state for each vertex.
		// The QP is generated such that the control input u=0 is always feasible.
		std::vector<DynamicVector> feasible_x(num_vertices(g));

		// Populate edge properties before vertex properties,
		// because making feasible state bounds relies on A and b.
		for (auto e : edges(g))
		{
			auto const& sz = qp.size();
			auto const u = source(e, g);
			auto const v = target(e, g);

			qp.A(e, rand_matrix.generate(sz.nx(v), sz.nx(u)));
			qp.B(e, rand_matrix.generate(sz.nx(v), sz.nu(u)));
			qp.b(e, rand_vector.generate(sz.nx(v)));
		}
		
		for (auto v : vert)
		{
			auto const& sz = qp.size();

			blaze::SymmetricMatrix<DynamicMatrix> H(sz.nx(v) + sz.nu(v));
			makePositiveDefinite(H);

			qp.Q(v, submatrix(H, 0, 0, sz.nx(v), sz.nx(v)));
			qp.q(v, rand_vector.generate(sz.nx(v)));
			qp.C(v, rand_matrix.generate(sz.nc(v), sz.nx(v)));		

			if (auto const parent_edge = g.parentEdge(v))
			{
				auto const parent_vertex = source(*parent_edge, g);				
				feasible_x[v] = qp.A(*parent_edge) * feasible_x[parent_vertex]
					+ qp.b(*parent_edge);
			}
			else
			{
				feasible_x[v] = blaze::ZeroVector<Real>(sz.nx(v));
			}

			qp.lx(v, feasible_x[v] - abs(rand_vector.generate(sz.nx(v))));
			qp.ux(v, feasible_x[v] + abs(rand_vector.generate(sz.nx(v))));

			qp.ld(v, qp.C(v) * feasible_x[v] - abs(rand_vector.generate(sz.nc(v))));
			qp.ud(v, qp.C(v) * feasible_x[v] + abs(rand_vector.generate(sz.nc(v))));

			if (out_degree(v, g) > 0)
			{
				qp.R(v, submatrix(H, sz.nx(v), sz.nx(v), sz.nu(v), sz.nu(v)));
				qp.S(v, submatrix(H, sz.nx(v), 0, sz.nu(v), sz.nx(v)));
				qp.r(v, rand_vector.generate(sz.nu(v)));
				qp.D(v, rand_matrix.generate(sz.nc(v), sz.nu(v)));
				qp.lu(v, -abs(rand_vector.generate(sz.nu(v))));
				qp.uu(v, abs(rand_vector.generate(sz.nu(v))));
			}

			// {
			// 	DynamicMatrix const Z = rand_matrix.generate(sz.ns(), sz.ns());
			// 	put(qp.Zl(ctrans(Z) * Z);
			// }

			// {
			// 	DynamicMatrix const Z = rand_matrix.generate(sz.ns(), sz.ns());
			// 	put(qp.Zu(ctrans(Z) * Z);
			// }
				
			// put(qp.zl(rand_vector.generate(sz.ns()));
			// put(qp.zu(rand_vector.generate(sz.ns()));
		}
	}
}