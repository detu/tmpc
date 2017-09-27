# Generates test data for ExplicitRungeKutta4 integrator test.
# The test system is a pendulum with friction.
# The dynamic equations and parameters are taken from ACADO pendulum example:
# https://github.com/acado/acado/blob/master/examples/integrator/pendulum.cpp 

import casadi as cs
import numpy as np
import os
import sys


def ensure_dir_exist(dirname):
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise
        

if __name__ == '__main__':
    name_c = sys.argv[1]

    phi   = cs.MX.sym('phi')     # the angle phi
    dphi  = cs.MX.sym('dphi')    # the first derivative of phi w.r.t time
    u     = cs.MX.sym('u')       # a force acting on the pendulum
    l     = 1.                   # the length of the pendulum
    m     = 1.0                  # the mass of the pendulum
    g     = 9.81                 # the gravitational constant
    alpha = 2.0                  # frictional constant
    ts    = 0.01                 # time step
    x     = cs.vertcat(phi, dphi)     # state vector
    
    z = cs.sin(phi)
    dx = cs.vertcat(dphi, -(m*g/l)*z - alpha*dphi + u/(m*l))
    q = cs.vertcat(dphi ** 2, phi)   # quadrature term
    r = cs.vertcat(cs.cos(phi), u) # residual term
    
    
    #------------------------------
    # Generate C code for ODE model
    #------------------------------
    name = 'pendulum_ode'
    ode = cs.Function(name, [x, u], 
                          [cs.densify(dx), cs.densify(cs.jacobian(dx, x)), cs.densify(cs.jacobian(dx, u))], 
                          ['x0', 'u0'], ['xdot', 'A', 'B'])
    
    gen = cs.CodeGenerator(name_c, {'mex' : False, 'with_header' : True, 'with_mem' : True})
    gen.add(ode)
    gen.generate()

    print("Generated code written to {0}".format(name_c))