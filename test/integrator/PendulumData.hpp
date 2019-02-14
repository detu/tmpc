#pragma once

#include "PendulumOde.hpp"

#include <iosfwd>
#include <vector>


namespace tmpc :: testing
{
    struct TestPoint
    {
        double t;
        PendulumOdeBase::StateVector xdot;
        PendulumOdeBase::StateStateMatrix Aode;
        PendulumOdeBase::StateInputMatrix Bode;
        PendulumOdeBase::QuadVector q;
        PendulumOdeBase::QuadStateMatrix qA_ode;
        PendulumOdeBase::QuadInputMatrix qB_ode;
        PendulumOdeBase::StateVector x0;
        PendulumOdeBase::InputVector u;
        PendulumOdeBase::StateVector xplus;
        PendulumOdeBase::QuadVector qf;
        PendulumOdeBase::StateStateMatrix A;
        PendulumOdeBase::StateInputMatrix B;
        PendulumOdeBase::QuadStateMatrix qA;
        PendulumOdeBase::QuadInputMatrix qB;
        PendulumOdeBase::ResVector r;
        PendulumOdeBase::ResStateMatrix rA_ode;
        PendulumOdeBase::ResInputMatrix rB_ode;
        double cf;
        PendulumOdeBase::StateVector cA;
        PendulumOdeBase::InputVector cB;
        PendulumOdeBase::StateStateMatrix cQ;
        PendulumOdeBase::InputInputMatrix cR;
        PendulumOdeBase::StateInputMatrix cS;
        double const timeStep = 1e-2;

        friend std::istream& operator>>(std::istream& is, TestPoint& p);
    };


    std::vector<TestPoint> loadPendulumData();
}