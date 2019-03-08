#include "PendulumData.hpp"


#include <fstream>


namespace tmpc :: testing
{
    template <typename MT, bool SO>
	inline std::istream& operator>>(std::istream& is, blaze::Matrix<MT, SO>& m)
	{
		for (size_t i = 0; i < rows(m); ++i)
			for (size_t j = 0; j < columns(m); ++j)
				is >> (~m)(i, j);

		return is;
	}

	
    template <typename VT, bool TF>
	inline std::istream& operator>>(std::istream& is, blaze::Vector<VT, TF>& v)
	{
		for (size_t i = 0; i < size(v); ++i)
			is >> (~v)[i];

		return is;
	}


    std::istream& operator>>(std::istream& is, TestPoint& p)
    {
        // keys = ['t', 'x0', 'u', 'xdot', 'A_ode', 'B_ode', 'q', 'qA_ode', 'qB_ode', 'x_plus', 'A', 'B', 'qf', 'qA', 'qB']
        return is >> p.t >> p.x0 >> p.u >> p.xdot  >> p.Aode   >> p.Bode
                                        >> p.q     >> p.qA_ode >> p.qB_ode
                                        >> p.xplus >> p.A      >> p.B
                                        >> p.qf    >> p.qA     >> p.qB
                                        >> p.r	   >> p.rA_ode >> p.rB_ode
                                        >> p.cf    >> p.cA     >> p.cB
                                        >> p.cQ    >> p.cR     >> p.cS;
    };


    std::vector<TestPoint> loadPendulumData()
    {
        std::ifstream data_file {std::string(TEST_DATA_PATH) + "/rk4/pendulum.txt"};
        return std::vector<TestPoint> {std::istream_iterator<TestPoint>(data_file), std::istream_iterator<TestPoint>()};
    }
}