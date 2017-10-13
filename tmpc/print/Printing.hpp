#pragma once

#include <ostream>

namespace tmpc
{
	template <typename Matrix>
	class MatlabPrintWrapper
	{
	public:
		MatlabPrintWrapper(Matrix const& val) : val_(val) {}

		void print(std::ostream& os) const
		{
			if (val_.rows() > 0 && val_.cols() > 0)
				os << "[..." << std::endl << val_ << "]";
			else
				os << "zeros(" << val_.rows() << ", " << val_.cols() << ")";
		}

	private:
		Matrix const& val_;
	};

	template <typename Matrix>
	inline MatlabPrintWrapper<Matrix> AsMatlab(Matrix const& val)
	{
		return MatlabPrintWrapper<Matrix>(val);
	}

	template <typename Matrix>
	inline std::ostream& operator<<(std::ostream& os, MatlabPrintWrapper<Matrix> const& w)
	{
		w.print(os);
		return os;
	}

	/**
	 * \brief Print a multistage QP as a MATLAB script.
	 *
	 * TODO: move all QPs to a qp namespace. Then move the printing function to the QP namespace.
	 * Then if a Print*** function is called for a QP class, it will be the function from the qp namespace.
	 */
	template <typename QP>
	inline void PrintMultistageQpMatlab(std::ostream& os, QP const& qp, const std::string& var_name)
	{
		using std::endl;

		unsigned k = 1;

		for (auto const& stage : qp)
		{
			std::string const qp_i = var_name + "(" + std::to_string(k) + ")";

			os << qp_i << ".Q = " << AsMatlab(stage.Q()) << ";" << endl << endl;
			os << qp_i << ".R = " << AsMatlab(stage.R()) << ";" << endl << endl;
			os << qp_i << ".S = " << AsMatlab(stage.S()) << ";" << endl << endl;
			os << qp_i << ".q = " << AsMatlab(stage.q()) << ";" << endl << endl;
			os << qp_i << ".r = " << AsMatlab(stage.r()) << ";" << endl << endl;

			os << qp_i << ".A = " << AsMatlab(stage.A()) << ";" << endl << endl;
			os << qp_i << ".B = " << AsMatlab(stage.B()) << ";" << endl << endl;
			os << qp_i << ".b = " << AsMatlab(stage.b()) << ";" << endl << endl;

			os << qp_i << ".C = " << AsMatlab(stage.C()) << ";" << endl << endl;
			os << qp_i << ".D = " << AsMatlab(stage.D()) << ";" << endl << endl;
			os << qp_i << ".lbd = " << AsMatlab(stage.lbd()) << ";" << endl << endl;
			os << qp_i << ".ubd = " << AsMatlab(stage.ubd()) << ";" << endl << endl;

			os << qp_i << ".lbx = " << AsMatlab(stage.lbx()) << ";" << endl << endl;
			os << qp_i << ".ubx = " << AsMatlab(stage.ubx()) << ";" << endl << endl;
			os << qp_i << ".lbu = " << AsMatlab(stage.lbu()) << ";" << endl << endl;
			os << qp_i << ".ubu = " << AsMatlab(stage.ubu()) << ";" << endl << endl;

			os << qp_i << ".Zl = " << AsMatlab(stage.Zl()) << ";" << endl << endl;
			os << qp_i << ".Zu = " << AsMatlab(stage.Zu()) << ";" << endl << endl;
			os << qp_i << ".zl = " << AsMatlab(stage.zl()) << ";" << endl << endl;
			os << qp_i << ".zu = " << AsMatlab(stage.zu()) << ";" << endl << endl;

			os << qp_i << ".idxs = [ ..." << endl;
			for (auto i: stage.idxs())
				os << i << ", ";
			os << "]" << endl << endl;
			
			++k;
		}
	}
}
