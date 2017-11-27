#pragma once

#include <tmpc/qp/OcpQpBase.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print an OCP QP stage in a human-readable format.
	 */
	template <typename Stage>
	inline std::ostream& operator<<(std::ostream& os, OcpQpBase<Stage> const& stage)
	{
		using std::endl;

		os << "Q = " << endl << stage.Q() << endl << endl;
		os << "R = " << endl << stage.R() << endl << endl;
		os << "S = " << endl << stage.S() << endl << endl;
		os << "q = " << endl << trans(stage.q()) << endl << endl;
		os << "r = " << endl << trans(stage.r()) << endl << endl;

		os << "A = " << endl << stage.A() << endl << endl;
		os << "B = " << endl << stage.B() << endl << endl;
		os << "b = " << endl << trans(stage.b()) << endl << endl;

		os << "C = " << endl << stage.C() << endl << endl;
		os << "D = " << endl << stage.D() << endl << endl;
		os << "lbd = " << endl << trans(stage.lbd()) << endl << endl;
		os << "ubd = " << endl << trans(stage.ubd()) << endl << endl;

		os << "lbx = " << endl << trans(stage.lbx()) << endl << endl;
		os << "ubx = " << endl << trans(stage.ubx()) << endl << endl;
		os << "lbu = " << endl << trans(stage.lbu()) << endl << endl;
		os << "ubu = " << endl << trans(stage.ubu()) << endl << endl;

		os << "Zl = " << endl << stage.Zl() << endl << endl;
		os << "Zu = " << endl << stage.Zu() << endl << endl;
		os << "zl = " << endl << trans(stage.zl()) << endl << endl;
		os << "zu = " << endl << trans(stage.zu()) << endl << endl;

		os << "idxs = " << endl;
		for (auto i: stage.idxs())
			os << i << " ";
		os << endl;

		return os;
	}
}
