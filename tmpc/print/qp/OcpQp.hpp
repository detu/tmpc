#pragma once

#include <tmpc/qp/OcpQp.hpp>

#include <ostream>


namespace tmpc
{
	/**
	 * \brief Print an OCP QP stage in a human-readable format.
	 */
	template <OcpQp Qp>
	inline std::ostream& operator<<(std::ostream& os, Qp const& qp)
	{
		using std::endl;

		for (auto v : vertices(qp.graph()))
		{
			bool const is_branch = out_degree(v, qp.graph()) > 0;

			os << "Q[" << v << "] = " << endl << qp.Q(v) << endl;
			if (is_branch) os << "R[" << v << "] = " << endl << qp.R(v) << endl;
			if (is_branch) os << "S[" << v << "] = " << endl << qp.S(v) << endl;
			os << "q[" << v << "] = " << endl << trans(qp.q(v)) << endl;
			if (is_branch) os << "r[" << v << "] = " << endl << trans(qp.r(v)) << endl;

			os << "C[" << v << "] = " << endl << qp.C(v) << endl;
			if (is_branch) os << "D[" << v << "] = " << endl << qp.D(v) << endl;
			os << "ld[" << v << "] = " << endl << trans(qp.ld(v)) << endl;
			os << "ud[" << v << "] = " << endl << trans(qp.ud(v)) << endl;

			os << "lx[" << v << "] = " << endl << trans(qp.lx(v)) << endl;
			os << "ux[" << v << "] = " << endl << trans(qp.ux(v)) << endl;
			if (is_branch) os << "lu[" << v << "] = " << endl << trans(qp.lu(v)) << endl;
			if (is_branch) os << "uu[" << v << "] = " << endl << trans(qp.uu(v)) << endl;

			// os << "Zl = " << endl << stage.Zl() << endl << endl;
			// os << "Zu = " << endl << stage.Zu() << endl << endl;
			// os << "zl = " << endl << trans(stage.zl()) << endl << endl;
			// os << "zu = " << endl << trans(stage.zu()) << endl << endl;

			// os << "idxs = " << endl;
			// for (auto i: stage.idxs())
			// 	os << i << " ";
			// os << endl;
		}

		for (auto e : edges(qp.graph()))
		{
			os << "A[" << e << "] = " << endl << qp.A(e) << endl;
			os << "B[" << e << "] = " << endl << qp.B(e) << endl;
			os << "b[" << e << "] = " << endl << trans(qp.b(e)) << endl;
		}

		return os;
	}
}
