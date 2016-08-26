# Generates c function for RK4 integrator.

import casadi as cs
import numpy as np
import os

t     = cs.MX.sym('t')       # time
x1    = cs.MX.sym('x1')      # state x1
x2    = cs.MX.sym('x2')      # state x2
u     = cs.MX.sym('u')       # control u
mu    = 0.5                  # abstract ode parameter
x     = cs.vertcat(x1, x2)   # state vector

f = cs.vertcat(x2 + u*(mu+(1-mu)*x1), x1 + u*(mu-4*(1-mu)*x2))

#------------------------------
# Generate C code for ODE model
#------------------------------
name = 'chen1998_ode'
ode_AB = cs.Function(name + "_AB", [t, x, u],
                      [cs.densify(f), cs.densify(cs.jacobian(f, x)), cs.densify(cs.jacobian(f, u))], ['t', 'x0', 'u0'], ['xdot', 'A', 'B'])
gen = cs.CodeGenerator({'mex' : False, 'with_header' : True})
gen.add(ode_AB)
name_c = '{0}_generated.c'.format(name)
name_h = '{0}_generated.h'.format(name)
gen.generate(name_c)
