/*
 * RK4.hpp
 *
 *  Created on: May 8, 2016
 *      Author: kotlyar
 */

#pragma once

#include <Eigen/Dense>

#include "../core/matrix.hpp"

namespace tmpc
{
	/**
	 * \defgroup integrators Integrators
	 */

	/**
	 * \brief Explicit Runge-Kutta 4 integrator.
	 * \ingroup integrators
	 */
	class RK4
	{
	public:
		RK4(double time_step) : _timeStep(time_step) {}

		// TODO: should it be an argument of the integrate() function instead?
		double timeStep() const noexcept { return _timeStep; }

	private:
		double const _timeStep;
	};

	template <typename ODE, typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix>
	void integrate(RK4 const& integrator, ODE const& ode, double t0, StateVector0_ const& x0, InputVector_ const& u, StateVector1_& x_next, AMatrix& A, BMatrix& B)
	{
		auto constexpr NX = rows<StateVector0_>();
		auto constexpr NU = rows<InputVector_ >();

		typedef Eigen::Matrix<double, NX,  1> StateVector;
		typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
		typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;

		StateVector k1, k2, k3, k4;
		StateStateMatrix A1, A2, A3, A4;
		StateInputMatrix B1, B2, B3, B4;
		auto const h = integrator.timeStep();

		// Calculating next state
		ode(t0,          x0                , u, k1, A1, B1);
		ode(t0 + h / 2., x0 + k1 * (h / 2.), u, k2, A2, B2);
		ode(t0 + h / 2., x0 + k2 * (h / 2.), u, k3, A3, B3);
		ode(t0 + h,      x0 + k3 * h       , u, k4, A4, B4);

		x_next = x0 + (k1 + 2. * k2 + 2. * k3 + k4) * (h / 6.);

		// Calculating sensitivities
		auto const& A1_bar =      A1;							auto const& B1_bar =      B1;
		auto const  A2_bar = eval(A2 + (h / 2.) * A2 * A1_bar);	auto const  B2_bar = eval(B2 + (h / 2.) * A2 * B1_bar);
		auto const  A3_bar = eval(A3 + (h / 2.) * A3 * A2_bar);	auto const  B3_bar = eval(B3 + (h / 2.) * A3 * B2_bar);
		auto const  A4_bar =      A4 +  h       * A4 * A3_bar ;	auto const  B4_bar =      B4 +  h       * A4 * B3_bar ;

		A = identity<StateStateMatrix>() + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					           (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);
	}

	template <typename ODE,
		typename InitialStateVector, typename InputVector,
		typename InitialStateSensitivityVector, typename InputSensitivityVector,
		typename FinalStateVector, typename FinalStateSensitivityVector>
	void integrate(RK4 const& integrator, ODE const& ode, double t0, InitialStateVector const& x0, InputVector const& u,
		InitialStateSensitivityVector const& seed_x0, InputSensitivityVector const& seed_u, FinalStateVector& xf, FinalStateSensitivityVector& sens_xf)
	{
		Eigen::Matrix<double, rows<InitialStateVector>(),  1> k1, k2, k3, k4, sens_k1, sens_k2, sens_k3, sens_k4;
		auto const h = integrator.timeStep();

		// Calculating next state
		ode(t0,          x0                , u, seed_x0                     , seed_u, k1, sens_k1);
		ode(t0 + h / 2., x0 + k1 * (h / 2.), u, seed_x0 + sens_k1 * (h / 2.), seed_u, k2, sens_k2);
		ode(t0 + h / 2., x0 + k2 * (h / 2.), u, seed_x0 + sens_k2 * (h / 2.), seed_u, k3, sens_k3);
		ode(t0 + h,      x0 + k3 *  h      , u, seed_x0 + sens_k3 *  h,       seed_u, k4, sens_k4);

		xf      =      x0 + (     k1 + 2. *      k2 + 2. *      k3 +      k4) * (h / 6.);
		sens_xf = seed_x0 + (sens_k1 + 2. * sens_k2 + 2. * sens_k3 + sens_k4) * (h / 6.);
	}

	//
	// Sensitivity-free version of the integrate() function.
	//
	template <typename ODE, typename StateVector, typename InputVector>
	decltype(auto) integrate(RK4 const& integrator, ODE const& ode, double t0, StateVector const& x0, InputVector const& u)
	{
		auto const h = integrator.timeStep();

		// Calculating next state
		auto const k1 = eval(ode(t0,          x0                , u));
		auto const k2 = eval(ode(t0 + h / 2., x0 + k1 * (h / 2.), u));
		auto const k3 = eval(ode(t0 + h / 2., x0 + k2 * (h / 2.), u));
		auto const k4 = eval(ode(t0 + h,      x0 + k3 * h       , u));

		return eval(x0 + (k1 + 2. * k2 + 2. * k3 + k4) * (h / 6.));
	}

	/**
	 * \brief Integration including a quadrature.
	 * \ingroup integrators
	 *
	 *
	 * Evaluates the following ODE using explicit Runge-Kutta 4 method:
	 * \f{eqnarray*}{
	 * \dot{\mathbf{x}} = f(t,\mathbf{x},\mathbf{u})\\
	 * \dot{\mathbf{q}} = h(t,\mathbf{x},\mathbf{u})
	 * \f}
	 * on time interval \f$[t_0,t_0+h]\f$ subject to initial conditions
	 * \f{align*}{
	 * \mathbf{x}(t_0) & =\mathbf{x}_0\\
	 * \mathbf{q}(t_0) & =0
	 * \f}
	 * \param[in] integrator -- integrator
	 * \param[in] ode -- ODE to integrate. ode(t, x, u, k, A, B, q, qA, qB) must be a valid expression (TODO: document ode input and output).
	 * \param[in] t0 -- start of the integration interval \f$t_0\f$
	 * \param[in] x0 -- state of the system at the beginning of integration interval \f$\mathbf{x}_0\f$
	 * \param[in] u -- input of the system which is assumed to be constant for the entire integration interval \f$\mathbf{u}\f$
	 * \param[out] x_next -- state of the system which at the end of the integration interval \f$\mathbf{x}(t_0+h)\f$
	 * \param[out] A -- sensitivity of final state w.r.t. initial state \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] B -- sensitivity of final state w.r.t. input \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] qf -- quadrature state at the end of the integration interval \f$\mathbf{q}(t_0+h)\f$
	 * \param[out] qA -- sensitivity of final quadrature state w.r.t. initial state \f$\frac{\mathrm{d}\mathbf{q}(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] qB -- sensitivity of final quadrature state w.r.t. input \f$\frac{\mathrm{d}\mathbf{q}(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 */
	template <typename ODE, typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix,
		typename QuadVector_, typename QuadStateMatrix_, typename QuadInputMatrix_>
	void integrate(RK4 const& integrator, ODE const& ode, double t0, StateVector0_ const& x0, InputVector_ const& u,
			StateVector1_& x_next, AMatrix& A, BMatrix& B, QuadVector_& qf, QuadStateMatrix_& qA, QuadInputMatrix_& qB)
	{
		auto constexpr NX = rows<StateVector0_>();
		auto constexpr NU = rows<InputVector_ >();
		auto constexpr NQ = rows<QuadVector_  >();

		typedef Eigen::Matrix<double, NX,  1> StateVector;
		typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
		typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;

		StateVector k1, k2, k3, k4;
		StateStateMatrix A1, A2, A3, A4;
		StateInputMatrix B1, B2, B3, B4;

		Eigen::Matrix<double, NQ,  1> dq1, dq2, dq3, dq4;
		Eigen::Matrix<double, NQ, NX> qA1, qA2, qA3, qA4;
		Eigen::Matrix<double, NQ, NU> qB1, qB2, qB3, qB4;
		auto const h = integrator.timeStep();

		// Calculating next state and quadrature
		ode(t0,          x0                , u, k1, A1, B1, dq1, qA1, qB1);
		ode(t0 + h / 2., x0 + k1 * (h / 2.), u, k2, A2, B2, dq2, qA2, qB2);
		ode(t0 + h / 2., x0 + k2 * (h / 2.), u, k3, A3, B3, dq3, qA3, qB3);
		ode(t0 + h,      x0 + k3 * h       , u, k4, A4, B4, dq4, qA4, qB4);

		x_next = x0 + ( k1 + 2. *  k2 + 2. *  k3 +  k4) * (h / 6.);
		qf     =      (dq1 + 2. * dq2 + 2. * dq3 + dq4) * (h / 6.);

		// Calculating sensitivities
		auto const& A1_bar =      A1;							auto const& B1_bar =      B1;
		auto const  A2_bar = eval(A2 + (h / 2.) * A2 * A1_bar);	auto const  B2_bar = eval(B2 + (h / 2.) * A2 * B1_bar);
		auto const  A3_bar = eval(A3 + (h / 2.) * A3 * A2_bar);	auto const  B3_bar = eval(B3 + (h / 2.) * A3 * B2_bar);
		auto const  A4_bar =      A4 +  h       * A4 * A3_bar ;	auto const  B4_bar =      B4 +  h       * A4 * B3_bar ;

		A = identity<StateStateMatrix>() + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					           (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);

		auto const& qA1_bar =      qA1;							    auto const& qB1_bar =      qB1;
		auto const  qA2_bar = eval(qA2 + (h / 2.) * qA2 * A1_bar);	auto const  qB2_bar = eval(qB2 + (h / 2.) * qA2 * B1_bar);
		auto const  qA3_bar = eval(qA3 + (h / 2.) * qA3 * A2_bar);	auto const  qB3_bar = eval(qB3 + (h / 2.) * qA3 * B2_bar);
		auto const  qA4_bar =      qA4 +  h       * qA4 * A3_bar ;	auto const  qB4_bar =      qB4 +  h       * qA4 * B3_bar ;

		qA = (h / 6.) * (qA1_bar + 2. * qA2_bar + 2. * qA3_bar + qA4_bar);
		qB = (h / 6.) * (qB1_bar + 2. * qB2_bar + 2. * qB3_bar + qB4_bar);
	}


	/*
	// Version of integrate() working with backward-mode ODE sensitivities. Not implemented properly yet.
	//
	template <typename Integrator, typename ODE,
		typename InitialStateVector, typename InputVector_, typename FinalStateVector, typename AMatrix, typename BMatrix>
	void integrate(Integrator const& integrator, ODE const& ode, double t0,
			InitialStateVector const& x0, InputVector_ const& u, FinalStateVector& xf, AMatrix& A, BMatrix& B)
	{
		auto constexpr NX = rows<InitialStateVector>();
		auto constexpr NU = rows<InputVector_      >();

		typedef Eigen::Matrix<double, NX,  1> StateVector;
		typedef Eigen::Matrix<double, NU,  1> InputVector;

		for (decltype(rows<InitialStateVector>()) i = 0; i < NX; ++i)
		{
			StateVector seed_x0 = zero<StateVector>();
			seed_x0(i) = 1.;

			auto col_i = col(A, i);
			integrator.Integrate(ode, t0, x0, u, seed_x0, zero<InputVector>(), xf, col_i);
		}

		for (decltype(rows<InputVector_      >()) i = 0; i < NU; ++i)
		{
			InputVector seed_u = zero<InputVector>();
			seed_u(i) = 1.;

			auto col_i = col(B, i);
			integrator.Integrate(ode, t0, x0, u, zero<StateVector>(), seed_u, xf, col_i);
		}
	}
	*/
}
