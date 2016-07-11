# Generates test data for RK4 integrator test.
# The test system is a pendulum with friction.
# The dynamic equations and parameters are taken from ACADO pendulum example:
# https://github.com/acado/acado/blob/master/examples/integrator/pendulum.cpp 

import casadi as cs
import matplotlib.pyplot as plt
import numpy as np
import os

t     = cs.MX.sym('t')       # time
phi   = cs.MX.sym('phi')     # the angle phi
dphi  = cs.MX.sym('dphi')    # the first derivative of phi w.r.t time
F     = cs.MX.sym('F')       # a force acting on the pendulum
l     = 1.                   # the length of the pendulum
m     = 1.0                  # the mass of the pendulum
g     = 9.81                 # the gravitational constant
alpha = 2.0                  # frictional constant
ts    = 0.01                 # time step
x     = cs.vertcat(phi, dphi)     # state vector

z = cs.sin(phi)
f = cs.vertcat(dphi, -(m*g/l)*z - alpha*dphi + F/(m*l))
ode = cs.Function('Pendulum_ODE', [x, F], [f])

integrator = cs.simpleRK(ode, 1)

x_plus = integrator(x0 = x, p = F, h = ts)['xf']
integrator_sens = cs.Function('Integrator_Sensitivities', [x, F], [x_plus, cs.jacobian(x_plus, x), cs.jacobian(x_plus, F)])

#------------------------------
# Generate C code for ODE model
#------------------------------
name = 'pendulum_ode'
#x_seed = cs.MX('x_seed', x.shape);
#u_seed = cs.MX('u_seed', F.shape);
ode_jac = cs.Function(name + "_jac", [t, cs.vertcat(x, F)], 
                      [cs.densify(f), cs.densify(cs.jacobian(f, cs.vertcat(x, F)))], ['t', 'z0'], ['xdot', 'J'])
ode_AB = cs.Function(name + "_AB", [t, x, F], 
                      [cs.densify(f), cs.densify(cs.jacobian(f, x)), cs.densify(cs.jacobian(f, F))], ['t', 'x0', 'u0'], ['xdot', 'A', 'B'])
#ode_sens = cs.Function(name + "_sens", [t, x, F, x_seed, u_seed], 
#                      [cs.densify(f), cs.densify(cs.jtimes(f, x, x_seed)), cs.densify(cs.jtimes(f, u, u_seed))], 
#                      ['t', 'x0', 'u0', 'x_seed', 'u_seed'], ['xdot', 'x_sens', 'u_sens'])
gen = cs.CodeGenerator({'mex' : False, 'with_header' : True})
gen.add(ode_jac)
gen.add(ode_AB)
name_c = '{0}_generated.c'.format(name)
name_h = '{0}_generated.h'.format(name)
gen.generate(name_c)
#os.rename(name_c, 'src/' + name_c)
#os.rename(name_h, 'src/' + name_h)

#------------------------------
# Generate test data
#------------------------------
x0 = cs.DM([1.0, 0.0])  # initial state

data = {'t' : [], 'z0' : [], 'xdot' : [], 'J_ode' : [], 'x_plus' : [], 'J' : []}
N = 600

for k in range(N):
    t_k = k * ts
    
    # input
    if t_k < 1.5:    
        u  = -10.0
    else:
        u = 0.0
                                        
    [x_plus, A, B] = integrator_sens(x0, u)
    
    z0  = cs.vertcat(x0, u)
    [xdot, J_ode] = ode_jac(t_k, z0)
    
    data['t'     ].append(t_k)
    data['z0'    ].append(z0)
    data['xdot'  ].append(xdot)
    data['J_ode' ].append(J_ode)
    data['x_plus'].append(x_plus)
    data['J'     ].append(cs.horzcat(A, B))
    x0 = x_plus
    
with open('data/rk4/pendulum.txt', 'w') as file:
    sep = ' '
    for k in range(N):
        np.array(data['t'     ][k]).tofile(file, sep)
        file.write('\t')
        np.array(data['z0'    ][k]).tofile(file, sep)
        file.write('\t')
        np.array(data['xdot'  ][k]).tofile(file, sep)
        file.write('\t')
        np.array(data['J_ode' ][k]).tofile(file, sep)
        file.write('\t')
        np.array(data['x_plus'][k]).tofile(file, sep)
        file.write('\t')
        np.array(data['J'     ][k]).tofile(file, sep)
        file.write('\n')
        
plt.subplot(2, 1, 1)
plt.step(cs.vertcat(data['t']), cs.transpose(cs.horzcat(*data['z0'])[ -1, :]))
plt.subplot(2, 1, 2)
plt.plot(cs.vertcat(data['t']), cs.transpose(cs.horzcat(*data['z0'])[: 2, :]), '.-')
plt.show()