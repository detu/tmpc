/*
 * RK4.hpp
 *
 *  Created on: May 8, 2016
 *      Author: kotlyar
 */

#pragma once

#include <tmpc/Matrix.hpp>

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
	void integrate(RK4 const& integrator, ODE const& ode, double t0,
			StateVector0_ const& x0,
			InputVector_ const& u,
			StateVector1_& x_next,
			AMatrix& A,
			BMatrix& B)
	{
		size_t constexpr NX = Size<StateVector0_>::value;
		size_t constexpr NU = Size<InputVector_ >::value;

		static_assert(Size<StateVector1_>::value == NX, "x_next must be of size NX");
		static_assert(Rows<AMatrix      >::value == NX && Columns<AMatrix      >::value == NX, "A must be of size NX*NX"    );
		static_assert(Rows<BMatrix      >::value == NX && Columns<BMatrix      >::value == NU, "B must be of size NX*NU"    );

		typedef StaticVector<double, NX, columnVector> StateVector;
		typedef StaticMatrix<double, NX, NX, columnMajor> StateStateMatrix;
		typedef StaticMatrix<double, NX, NU, columnMajor> StateInputMatrix;

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

		A = IdentityMatrix<double>(NX) + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					                 (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);
	}

	template <typename ODE,
		typename InitialStateVector, typename InputVector,
		typename InitialStateSensitivityVector, typename InputSensitivityVector,
		typename FinalStateVector, typename FinalStateSensitivityVector>
	void integrate(RK4 const& integrator, ODE const& ode, double t0, InitialStateVector const& x0, InputVector const& u,
		InitialStateSensitivityVector const& seed_x0, InputSensitivityVector const& seed_u, FinalStateVector& xf, FinalStateSensitivityVector& sens_xf)
	{
		StaticVector<double, Rows<InitialStateVector>::value> k1, k2, k3, k4, sens_k1, sens_k2, sens_k3, sens_k4;
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
	 * \param[in] ode -- ODE to integrate. ode(t, x, u, k, A, B, h, hA, hB) must be a valid expression (TODO: document ode input and output).
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
		auto constexpr NX = Rows<StateVector0_>::value;
		auto constexpr NU = Rows<InputVector_ >::value;
		auto constexpr NQ = Rows<QuadVector_  >::value;

		static_assert(Rows<StateVector1_   >::value == NX && Columns<StateVector1_   >::value ==  1, "x_next must be of size NX*1");
		static_assert(Rows<AMatrix         >::value == NX && Columns<AMatrix         >::value == NX, "A must be of size NX*NX"    );
		static_assert(Rows<BMatrix         >::value == NX && Columns<BMatrix         >::value == NU, "B must be of size NX*NU"    );
		static_assert(Rows<QuadVector_     >::value == NQ && Columns<QuadVector_     >::value ==  1, "qf must be of size NQ*1"    );
		static_assert(Rows<QuadStateMatrix_>::value == NQ && Columns<QuadStateMatrix_>::value == NX, "qA must be of size NQ*NX"   );
		static_assert(Rows<QuadInputMatrix_>::value == NQ && Columns<QuadInputMatrix_>::value == NU, "qB must be of size NQ*NU"   );

		typedef StaticVector<double, NX> StateVector;
		typedef StaticMatrix<double, NX, NX, columnMajor> StateStateMatrix;
		typedef StaticMatrix<double, NX, NU, columnMajor> StateInputMatrix;

		StateVector k1, k2, k3, k4;
		StateStateMatrix A1, A2, A3, A4;
		StateInputMatrix B1, B2, B3, B4;

		StaticVector<double, NQ> dq1, dq2, dq3, dq4;
		StaticMatrix<double, NQ, NX, columnMajor> qA1, qA2, qA3, qA4;
		StaticMatrix<double, NQ, NU, columnMajor> qB1, qB2, qB3, qB4;
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

		A = IdentityMatrix<double>(NX) + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					                 (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);

		auto const& qA1_bar =      qA1;							    auto const& qB1_bar =      qB1;
		auto const  qA2_bar = eval(qA2 + (h / 2.) * qA2 * A1_bar);	auto const  qB2_bar = eval(qB2 + (h / 2.) * qA2 * B1_bar);
		auto const  qA3_bar = eval(qA3 + (h / 2.) * qA3 * A2_bar);	auto const  qB3_bar = eval(qB3 + (h / 2.) * qA3 * B2_bar);
		auto const  qA4_bar =      qA4 +  h       * qA4 * A3_bar ;	auto const  qB4_bar =      qB4 +  h       * qA4 * B3_bar ;

		qA = (h / 6.) * (qA1_bar + 2. * qA2_bar + 2. * qA3_bar + qA4_bar);
		qB = (h / 6.) * (qB1_bar + 2. * qB2_bar + 2. * qB3_bar + qB4_bar);
	}

	/**
	 * \brief Integration including a cost term.
	 * \ingroup integrators
	 *
	 *
	 * Evaluates the following ODE using explicit Runge-Kutta 4 method:
	 * \f{align*}{
	 * \dot{\mathbf{x}} &= f(t,\mathbf{x},\mathbf{u})\\
	 * \dot{c} &= \frac{1}{2} \Vert r(t,\mathbf{x},\mathbf{u}) \Vert^2
	 * \f}
	 * on time interval \f$[t_0,t_0+h]\f$ subject to initial conditions
	 * \f{align*}{
	 * \mathbf{x}(t_0) & =\mathbf{x}_0\\
	 * c(t_0) & =0
	 * \f}
	 * \param[in] integrator -- integrator
	 * \param[in] ode -- ODE to integrate. ode(t, x, u, k, A, B, r, rA, rB) must be a valid expression (TODO: document ode input and output).
	 * \param[in] t0 -- start of the integration interval \f$t_0\f$
	 * \param[in] x0 -- state of the system at the beginning of integration interval \f$\mathbf{x}_0\f$
	 * \param[in] u -- input of the system which is assumed to be constant for the entire integration interval \f$\mathbf{u}\f$
	 * \param[out] x_next -- state of the system which at the end of the integration interval \f$\mathbf{x}(t_0+h)\f$
	 * \param[out] A -- sensitivity of final state w.r.t. initial state \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] B -- sensitivity of final state w.r.t. input \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] cf -- cost value at the end of the integration interval \f$c(t_0+h)\f$
	 * \param[out] cA -- sensitivity of final cost w.r.t. initial state \f$\frac{\mathrm{d}c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] cB -- sensitivity of final cost w.r.t. input \f$\frac{\mathrm{d}c(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] cQ -- second derivative of final cost w.r.t. initial state \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)^2}\f$
	 * \param[out] cR -- second derivative of final cost w.r.t. input \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{u}^2}\f$
	 * \param[out] cS -- second derivative of final cost w.r.t. initial state and input \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)\mathrm{d}\mathbf{u}}\f$
	 */
	template <typename ODE, typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix,
		typename StateVector2_,	typename InputVector2_,	typename QMatrix, typename RMatrix, typename SMatrix>
	void integrate(RK4 const& integrator, ODE const& ode, double t0, StateVector0_ const& x0, InputVector_ const& u,
			StateVector1_& x_next, AMatrix& A, BMatrix& B,
			double& cf, StateVector2_& cA, InputVector2_& cB, QMatrix& cQ, RMatrix& cR, SMatrix& cS)
	{
		auto constexpr NX = Rows<StateVector0_>::value;
		auto constexpr NU = Rows<InputVector_ >::value;
		auto constexpr NR = ODE::NR;

		static_assert(Rows<StateVector1_>::value == NX && Columns<StateVector1_>::value ==  1, "x_next must be of size NX*1");
		static_assert(Rows<AMatrix      >::value == NX && Columns<AMatrix      >::value == NX, "A must be of size NX*NX"    );
		static_assert(Rows<BMatrix      >::value == NX && Columns<BMatrix      >::value == NU, "B must be of size NX*NU"    );
		static_assert(Rows<StateVector2_>::value == NX && Columns<StateVector2_>::value ==  1, "cA must be of size NX*1"    );
		static_assert(Rows<InputVector2_>::value == NU && Columns<StateVector2_>::value ==  1, "cB must be of size NU*1"    );
		static_assert(Rows<QMatrix      >::value == NX && Columns<QMatrix      >::value == NX, "cQ must be of size NX*NX"   );
		static_assert(Rows<RMatrix      >::value == NU && Columns<RMatrix      >::value == NU, "cR must be of size NX*NX"   );
		static_assert(Rows<SMatrix      >::value == NX && Columns<SMatrix   >::value == NU, "cS must be of size NX*NU"   );

		typedef StaticVector<double, NX> StateVector;
		typedef StaticMatrix<double, NX, NX, columnMajor> StateStateMatrix;
		typedef StaticMatrix<double, NX, NU, columnMajor> StateInputMatrix;

		StateVector k1, k2, k3, k4;
		StateStateMatrix A1, A2, A3, A4;
		StateInputMatrix B1, B2, B3, B4;

		StaticVector<double, NR>  r1,  r2,  r3,  r4;
		StaticMatrix<double, NR, NX, columnMajor> rA1, rA2, rA3, rA4;
		StaticMatrix<double, NR, NU, columnMajor> rB1, rB2, rB3, rB4;
		auto const h = integrator.timeStep();

		// Calculating next state, quadrature and cost
		ode(t0,          x0                , u, k1, A1, B1, r1, rA1, rB1);
		ode(t0 + h / 2., x0 + k1 * (h / 2.), u, k2, A2, B2, r2, rA2, rB2);
		ode(t0 + h / 2., x0 + k2 * (h / 2.), u, k3, A3, B3, r3, rA3, rB3);
		ode(t0 + h,      x0 + k3 * h       , u, k4, A4, B4, r4, rA4, rB4);

		x_next =  x0 + (            k1  + 2. *              k2 + 2. *              k3 +              k4) * (h / 6.);
		cf     = 0.5 * (squaredNorm(r1) + 2. * squaredNorm(r2) + 2. * squaredNorm(r3) + squaredNorm(r4)) * (h / 6.);

		// Calculating sensitivities
		auto const& A1_bar =      A1;							auto const& B1_bar =      B1;
		auto const  A2_bar = eval(A2 + (h / 2.) * A2 * A1_bar);	auto const  B2_bar = eval(B2 + (h / 2.) * A2 * B1_bar);
		auto const  A3_bar = eval(A3 + (h / 2.) * A3 * A2_bar);	auto const  B3_bar = eval(B3 + (h / 2.) * A3 * B2_bar);
		auto const  A4_bar =      A4 +  h       * A4 * A3_bar ;	auto const  B4_bar =      B4 +  h       * A4 * B3_bar ;

		A = IdentityMatrix<double>(NX) + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					                 (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);

		auto const& rA1_bar =      rA1;							    auto const& rB1_bar =      rB1;
		auto const  rA2_bar = eval(rA2 + (h / 2.) * rA2 * A1_bar);	auto const  rB2_bar = eval(rB2 + (h / 2.) * rA2 * B1_bar);
		auto const  rA3_bar = eval(rA3 + (h / 2.) * rA3 * A2_bar);	auto const  rB3_bar = eval(rB3 + (h / 2.) * rA3 * B2_bar);
		auto const  rA4_bar = eval(rA4 +  h       * rA4 * A3_bar);	auto const  rB4_bar = eval(rB4 +  h       * rA4 * B3_bar);

		cA = (h / 6.) * (trans(r1) * rA1_bar + 2. * trans(r2) * rA2_bar + 2. * trans(r3) * rA3_bar + trans(r4) * rA4_bar);
		cB = (h / 6.) * (trans(r1) * rB1_bar + 2. * trans(r2) * rB2_bar + 2. * trans(r3) * rB3_bar + trans(r4) * rB4_bar);

		// Gauss-Newton approximation of the Hessian.
		cQ = (h / 6.) * (trans(rA1_bar) * rA1_bar + 2. * trans(rA2_bar) * rA2_bar + 2. * trans(rA3_bar) * rA3_bar + trans(rA4_bar) * rA4_bar);
		cR = (h / 6.) * (trans(rB1_bar) * rB1_bar + 2. * trans(rB2_bar) * rB2_bar + 2. * trans(rB3_bar) * rB3_bar + trans(rB4_bar) * rB4_bar);
		cS = (h / 6.) * (trans(rA1_bar) * rB1_bar + 2. * trans(rA2_bar) * rB2_bar + 2. * trans(rA3_bar) * rB3_bar + trans(rA4_bar) * rB4_bar);
	}

	/**
	 * \brief Integration including a quadrature and a cost term.
	 * \ingroup integrators
	 *
	 *
	 * Evaluates the following ODE using explicit Runge-Kutta 4 method:
	 * \f{align*}{
	 * \dot{\mathbf{x}} &= f(t,\mathbf{x},\mathbf{u})\\
	 * \dot{\mathbf{q}} &= h(t,\mathbf{x},\mathbf{u})\\
	 * \dot{c} &= \frac{1}{2} \Vert r(t,\mathbf{x},\mathbf{u}) \Vert^2
	 * \f}
	 * on time interval \f$[t_0,t_0+h]\f$ subject to initial conditions
	 * \f{align*}{
	 * \mathbf{x}(t_0) & =\mathbf{x}_0\\
	 * \mathbf{q}(t_0) & =0\\
	 * c(t_0) & =0
	 * \f}
	 * \param[in] integrator -- integrator
	 * \param[in] ode -- ODE to integrate. ode(t, x, u, k, A, B, h, hA, hB, r, rA, rB) must be a valid expression (TODO: document ode input and output).
	 * \param[in] t0 -- start of the integration interval \f$t_0\f$
	 * \param[in] x0 -- state of the system at the beginning of integration interval \f$\mathbf{x}_0\f$
	 * \param[in] u -- input of the system which is assumed to be constant for the entire integration interval \f$\mathbf{u}\f$
	 * \param[out] x_next -- state of the system which at the end of the integration interval \f$\mathbf{x}(t_0+h)\f$
	 * \param[out] A -- sensitivity of final state w.r.t. initial state \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] B -- sensitivity of final state w.r.t. input \f$\frac{\mathrm{d}\mathbf{x}(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] qf -- quadrature state at the end of the integration interval \f$\mathbf{q}(t_0+h)\f$
	 * \param[out] qA -- sensitivity of final quadrature state w.r.t. initial state \f$\frac{\mathrm{d}\mathbf{q}(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] qB -- sensitivity of final quadrature state w.r.t. input \f$\frac{\mathrm{d}\mathbf{q}(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] cf -- cost value at the end of the integration interval \f$c(t_0+h)\f$
	 * \param[out] cA -- sensitivity of final cost w.r.t. initial state \f$\frac{\mathrm{d}c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)}\f$
	 * \param[out] cB -- sensitivity of final cost w.r.t. input \f$\frac{\mathrm{d}c(t_0+h)}{\mathrm{d}\mathbf{u}}\f$
	 * \param[out] cQ -- second derivative of final cost w.r.t. initial state \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)^2}\f$
	 * \param[out] cR -- second derivative of final cost w.r.t. input \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{u}^2}\f$
	 * \param[out] cS -- second derivative of final cost w.r.t. initial state and input \f$\frac{\mathrm{d}^2 c(t_0+h)}{\mathrm{d}\mathbf{x}(t_0)\mathrm{d}\mathbf{u}}\f$
	 */
	template <typename ODE, typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix,
		typename QuadVector_, typename QuadStateMatrix_, typename QuadInputMatrix_,
		typename StateVector2_,	typename InputVector2_,	typename QMatrix, typename RMatrix, typename SMatrix>
	void integrate(RK4 const& integrator, ODE const& ode, double t0,
			StateVector0_ const& x0,
			InputVector_ const& u,
			StateVector1_& x_next,
			AMatrix& A,
			BMatrix& B,
			QuadVector_& qf,
			QuadStateMatrix_& qA,
			QuadInputMatrix_& qB,
			double& cf,
			StateVector2_& cA,
			InputVector2_& cB,
			QMatrix& cQ,
			RMatrix& cR,
			SMatrix& cS)
	{
		auto constexpr NX = Rows<StateVector0_>::value;
		auto constexpr NU = Rows<InputVector_ >::value;
		auto constexpr NQ = Rows<QuadVector_  >::value;
		auto constexpr NR = ODE::NR;

		static_assert(Rows<StateVector1_   >::value == NX && Columns<StateVector1_   >::value ==  1, "x_next must be of size NX*1");
		static_assert(Rows<AMatrix         >::value == NX && Columns<AMatrix         >::value == NX, "A must be of size NX*NX"    );
		static_assert(Rows<BMatrix         >::value == NX && Columns<BMatrix         >::value == NU, "B must be of size NX*NU"    );
		static_assert(Rows<QuadVector_     >::value == NQ && Columns<QuadVector_     >::value ==  1, "qf must be of size NQ*1"    );
		static_assert(Rows<QuadStateMatrix_>::value == NQ && Columns<QuadStateMatrix_>::value == NX, "qA must be of size NQ*NX"   );
		static_assert(Rows<QuadInputMatrix_>::value == NQ && Columns<QuadInputMatrix_>::value == NU, "qB must be of size NQ*NU"   );
		static_assert(Rows<StateVector2_   >::value == NX && Columns<StateVector2_   >::value ==  1, "cA must be of size NX*1"    );
		static_assert(Rows<InputVector2_   >::value == NU && Columns<StateVector2_   >::value ==  1, "cB must be of size NU*1"    );
		static_assert(Rows<QMatrix         >::value == NX && Columns<QMatrix         >::value == NX, "cQ must be of size NX*NX"   );
		static_assert(Rows<RMatrix         >::value == NU && Columns<RMatrix         >::value == NU, "cR must be of size NX*NX"   );
		static_assert(Rows<SMatrix         >::value == NX && Columns<SMatrix         >::value == NU, "cS must be of size NX*NU"   );

		typedef StaticVector<double, NX> StateVector;
		typedef StaticMatrix<double, NX, NX, columnMajor> StateStateMatrix;
		typedef StaticMatrix<double, NX, NU, columnMajor> StateInputMatrix;

		StateVector k1, k2, k3, k4;
		StateStateMatrix A1, A2, A3, A4;
		StateInputMatrix B1, B2, B3, B4;

		StaticVector<double, NQ> dq1, dq2, dq3, dq4;
		StaticMatrix<double, NQ, NX, columnMajor> qA1, qA2, qA3, qA4;
		StaticMatrix<double, NQ, NU, columnMajor> qB1, qB2, qB3, qB4;
		StaticVector<double, NR>  r1,  r2,  r3,  r4;
		StaticMatrix<double, NR, NX, columnMajor> rA1, rA2, rA3, rA4;
		StaticMatrix<double, NR, NU, columnMajor> rB1, rB2, rB3, rB4;
		auto const h = integrator.timeStep();

		// Calculating next state, quadrature and cost
		ode(t0,          x0                , u, k1, A1, B1, dq1, qA1, qB1, r1, rA1, rB1);
		ode(t0 + h / 2., x0 + k1 * (h / 2.), u, k2, A2, B2, dq2, qA2, qB2, r2, rA2, rB2);
		ode(t0 + h / 2., x0 + k2 * (h / 2.), u, k3, A3, B3, dq3, qA3, qB3, r3, rA3, rB3);
		ode(t0 + h,      x0 + k3 * h       , u, k4, A4, B4, dq4, qA4, qB4, r4, rA4, rB4);

		x_next =  x0 + (             k1  + 2. *               k2 + 2. *               k3 +               k4) * (h / 6.);
		qf     =       (             dq1 + 2. *              dq2 + 2. *              dq3 +              dq4) * (h / 6.);
		cf     = 0.5 * (squaredNorm(r1) + 2. * squaredNorm(r2) + 2. * squaredNorm(r3) + squaredNorm(r4)) * (h / 6.);

		// Calculating sensitivities
		auto const& A1_bar =      A1;							auto const& B1_bar =      B1;
		auto const  A2_bar = eval(A2 + (h / 2.) * A2 * A1_bar);	auto const  B2_bar = eval(B2 + (h / 2.) * A2 * B1_bar);
		auto const  A3_bar = eval(A3 + (h / 2.) * A3 * A2_bar);	auto const  B3_bar = eval(B3 + (h / 2.) * A3 * B2_bar);
		auto const  A4_bar =      A4 +  h       * A4 * A3_bar ;	auto const  B4_bar =      B4 +  h       * A4 * B3_bar ;

		A = IdentityMatrix<double>(NX) + (h / 6.) * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
		B = 					                 (h / 6.) * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);

		auto const& qA1_bar =      qA1;							    auto const& qB1_bar =      qB1;
		auto const  qA2_bar = eval(qA2 + (h / 2.) * qA2 * A1_bar);	auto const  qB2_bar = eval(qB2 + (h / 2.) * qA2 * B1_bar);
		auto const  qA3_bar = eval(qA3 + (h / 2.) * qA3 * A2_bar);	auto const  qB3_bar = eval(qB3 + (h / 2.) * qA3 * B2_bar);
		auto const  qA4_bar =      qA4 +  h       * qA4 * A3_bar ;	auto const  qB4_bar =      qB4 +  h       * qA4 * B3_bar ;

		qA = (h / 6.) * (qA1_bar + 2. * qA2_bar + 2. * qA3_bar + qA4_bar);
		qB = (h / 6.) * (qB1_bar + 2. * qB2_bar + 2. * qB3_bar + qB4_bar);

		auto const& rA1_bar =      rA1;							    auto const& rB1_bar =      rB1;
		auto const  rA2_bar = eval(rA2 + (h / 2.) * rA2 * A1_bar);	auto const  rB2_bar = eval(rB2 + (h / 2.) * rA2 * B1_bar);
		auto const  rA3_bar = eval(rA3 + (h / 2.) * rA3 * A2_bar);	auto const  rB3_bar = eval(rB3 + (h / 2.) * rA3 * B2_bar);
		auto const  rA4_bar = eval(rA4 +  h       * rA4 * A3_bar);	auto const  rB4_bar = eval(rB4 +  h       * rA4 * B3_bar);

		cA = (h / 6.) * (trans(r1) * rA1_bar + 2. * trans(r2) * rA2_bar + 2. * trans(r3) * rA3_bar + trans(r4) * rA4_bar);
		cB = (h / 6.) * (trans(r1) * rB1_bar + 2. * trans(r2) * rB2_bar + 2. * trans(r3) * rB3_bar + trans(r4) * rB4_bar);

		// Gauss-Newton approximation of the Hessian.
		cQ = (h / 6.) * (trans(rA1_bar) * rA1_bar + 2. * trans(rA2_bar) * rA2_bar + 2. * trans(rA3_bar) * rA3_bar + trans(rA4_bar) * rA4_bar);
		cR = (h / 6.) * (trans(rB1_bar) * rB1_bar + 2. * trans(rB2_bar) * rB2_bar + 2. * trans(rB3_bar) * rB3_bar + trans(rB4_bar) * rB4_bar);
		cS = (h / 6.) * (trans(rA1_bar) * rB1_bar + 2. * trans(rA2_bar) * rB2_bar + 2. * trans(rA3_bar) * rB3_bar + trans(rA4_bar) * rB4_bar);
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
