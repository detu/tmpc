#include "../include/core/mhe.hpp"
#include "../include/qp/CondensingSolver.hpp"
#include "../include/qp/HPMPCSolver.hpp"
#include "../include/integrator/rk4.hpp"

#include <Eigen/Dense>

#include <tmpc/Testing.hpp>

#include <type_traits>

struct ODE
{
	static unsigned const NX = 2;
	static unsigned const NU = 1;
	static unsigned const NW = 1;
	static unsigned const NY = 1;
	static unsigned const ND = 0;
	static unsigned const NDT = 0;

	typedef Eigen::Matrix<double, NU,  1> InputVector;
	typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
	typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;
	typedef Eigen::Matrix<double, NX, NW> StateDisturbanceMatrix;
	typedef Eigen::Matrix<double, NY, NX> OutputStateMatrix;
	typedef Eigen::Matrix<double, NY, NU> OutputInputMatrix;

	StateStateMatrix  const& A() const { return A_; }
	StateInputMatrix  const& B() const { return B_; }
	OutputStateMatrix const& C() const { return C_; }
	OutputInputMatrix const& D() const { return D_; }
	InputVector       const& u() const { return u_; }

	ODE(InputVector const& u)
	:	u_(u)
	{
		A_ << 1.,  1.,
			  0.,  1.;

		B_ << 0.5,
			  1. ;

		C_ << 1.,
			  0.;

		D_ << 0.;
	}

private:
	InputVector const& u_;
	StateStateMatrix A_;
	StateInputMatrix B_;
	OutputStateMatrix C_;
	OutputInputMatrix D_;
};

/**
 * \brief Sample trajectory estimation problem.
 *
 * The model is a double integrator with disturbance in input and output:
 * x^{+} = [1, 1; 0, 1] x + [0.5; 1] (u + w)
 * y     =       [1, 0] x +      [0]  u +     v
 */
class SampleTrajectoryEstimation
{
public:
	static unsigned const NX  = ODE::NX;
	static unsigned const NU  = ODE::NU;
	static unsigned const NW  = ODE::NW;
	static unsigned const NY  = ODE::NY;
	static unsigned const NC  = 0;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;

	SampleTrajectoryEstimation()
	{
		_x_min << -1., -1.;
		_x_max <<  1.,  1.;
		_u_min << -1.;
		_u_max <<  1.;
		_x_terminal_min << -1., -1.;
		_x_terminal_max <<  1.,  1.;

		x0_ << 0., 0.;
		invSigmaX0_ << 100.,   0.,
				         0., 100.;
		invSigmaW_  << 100.;
		invSigmaY_  << 100.;
	}

	ODE getODE(InputVector const& u) const
	{
		return ODE(u);
	}

	/**
	 * \brief Calculates stage cost of the MHE problem, its gradient and Hessian.
	 */
	template <
		typename StateVector, typename InputVector, typename DisturbanceVector, typename OutputVector,
		typename QMatrix, typename RMatrix, typename SMatrix, typename StateGradientVector, typename DisturbanceGradientVector>
	void StageCost(unsigned i,
		StateVector const& x, InputVector const& u, DisturbanceVector const& w, OutputVector const& y,
		QMatrix& Q, RMatrix& R, SMatrix& S,	StateGradientVector& q, DisturbanceGradientVector& r) const
	{
		Q = trans(ode_.C()) * invSigmaY_ * ode_.C();
		R = invSigmaW_;
		S = trans(ode_.C()) * invSigmaY_ * ode_.D();
		q = trans(ode_.C()) * invSigmaY_ * (ode_.C() * x + ode_.D() * u - y);
		r = invSigmaW_ * w;
	}

	StateVector const& getStateMin() const { return _x_min; }
	StateVector const& getStateMax() const { return _x_max; }

	InputVector const& getInputMin() const { return _u_min;	}
	InputVector const& getInputMax() const { return _u_max;	}

	/**
	 * \brief Calculates arrival cost of the MHE problem, its gradient and Hessian.
	 */
	template <typename StateVector, typename GradientVector, typename HessianMatrix>
	void ArrivalCost(const StateVector& x, GradientVector& g, HessianMatrix& H) const
	{
		H = invSigmaX0_;
		g = invSigmaX0_ * (x - x0_);
	}

	template <typename StateInputVector, typename ConstraintJacobianMatrix, typename ConstraintVector>
	void PathConstraints(unsigned i, const StateInputVector& z,
		ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{
	}

	template <typename StateVector, typename TerminalConstraintJacobianMatrix, typename TerminalConstraintVector>
	void TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
		TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const
	{
	}

private:
	StateVector _x_min;
	StateVector _x_max;
	StateVector _x_terminal_min;
	StateVector _x_terminal_max;
	InputVector _u_min;
	InputVector _u_max;

	StateVector x0_;	// Mean of the initial state estimate
	Eigen::Matrix<double, NX, NX> invSigmaX0_; // Inverse covariance of the initial state estimate
	Eigen::Matrix<double, NW, NW> invSigmaW_ ; // Inverse covariance of the disturbance
	Eigen::Matrix<double, NY, NY> invSigmaY_ ; // Inverse covariance of the measurement error
};

typedef SampleTrajectoryEstimation OCP;

class DiscreteTimeModel
{
public:
	double timeStep() const { return 1.; }
};

template <typename StateVector, typename InputVector, typename NextStateVector, typename AMatrix, typename BMatrix>
void integrate(DiscreteTimeModel const& integrator, ODE const& model, double t0, StateVector const& x0, InputVector const& u, NextStateVector& x_next,
		Eigen::MatrixBase<AMatrix>& A, Eigen::MatrixBase<BMatrix>& B)
{
	A = model.A();
	B = model.B();

	x_next = model.A() * x0 + model.B() * u;
}

template <typename RealtimeIteration>
class RealtimeIterationTest : public Test
{
public:
	RealtimeIterationTest(unsigned Nt = 2)
	:	_ocp(Nt)
	,	_qpSolver(Nt)
	,	_rti(_ocp, _integrator, _qpSolver, WorkingPoint(Nt, OCP::StateVector::Zero(), OCP::InputVector::Zero()))
	{
	}

protected:
	typedef typename RealtimeIteration::Integrator Integrator;
	//typedef tmpc::RealtimeIteration<OCP, Integrator, QPSolver> RealtimeIteration;
	typedef typename RealtimeIteration::WorkingPoint WorkingPoint;
	typedef typename RealtimeIteration::QPSolver QPSolver;

	OCP _ocp;
	Integrator _integrator;
	QPSolver _qpSolver;
	RealtimeIteration _rti;

	void Preparation()
	{
		_rti.Preparation();
	}

	OCP::InputVector Feedback(OCP::StateVector const& x0)
	{
		return _rti.Feedback(x0);
	}
};

typedef ::testing::Types<
		tmpc::RealtimeIteration<OCP, DiscreteTimeModel, tmpc::CondensingSolver>
,		tmpc::RealtimeIteration<OCP, DiscreteTimeModel, tmpc::HPMPCSolver     >
	> RTITypes;

TYPED_TEST_CASE(RealtimeIterationTest, RTITypes);

TYPED_TEST(RealtimeIterationTest, GivesCorrectU0)
{
	OCP::StateVector x;
	OCP::InputVector u;

	// Step 0
	{
		x << 1, 0;
		u = this->Feedback(x);
		this->Preparation();

		OCP::InputVector u_expected;
		u_expected << -0.690877362606266;

		EXPECT_TRUE(u.isApprox(u_expected, 1e-6));
	}

	// Step 1
	{
		{
			OCP::StateVector x;
			OCP::InputVector u;

			x << 0.654561318696867,	 -0.690877362606266;	u << 0.215679569867116;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(0).isApprox(x, 1e-6));
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_u(0).isApprox(u, 1e-6));

			x << 0.0715237410241597, -0.475197792739149;	u << 0.215679569867116;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(1).isApprox(x, 1e-6));
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_u(1).isApprox(u, 1e-6));

			x << 0.0715237410241597, -0.475197792739149;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(2).isApprox(x, 1e-6));
		}

		x << 0.654561318696867,	-0.690877362606266;
		u = this->Feedback(x);

		OCP::InputVector u_expected;
		u_expected << 0.218183;
		EXPECT_TRUE(u.isApprox(u_expected, 1e-5));
	}
}
