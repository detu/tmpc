# Generates test data for RK4 integrator test.
# The test system is a pendulum with friction.
# The dynamic equations and parameters are taken from ACADO pendulum example:
# https://github.com/acado/acado/blob/master/examples/integrator/pendulum.cpp 

import casadi as cs
import matplotlib.pyplot as plt
import numpy as np
import os

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
        
# Gauss-Newton RK4 integrator
def rk4gn(name, ode, options):
    t = cs.MX.sym('t')
    h = options['tf']
    f = cs.Function('f', [ode['x'], ode['p'], t], [ode['ode'], ode['quad'], ode['r']],
                    ['x', 'p', 't'], ['ode', 'quad', 'r'])
    
    x0 = cs.MX.sym('x0', ode['x'].shape)
    t0 = cs.MX.sym('t0')
    
    k1 = f(t = t0         , x = x0                       , p = u)
    k2 = f(t = t0 + h / 2., x = x0 + k1['ode'] * (h / 2.), p = u)
    k3 = f(t = t0 + h / 2., x = x0 + k2['ode'] * (h / 2.), p = u)
    k4 = f(t = t0 + h,      x = x0 + k3['ode'] *  h      , p = u)

    xf = x0 + (k1['ode' ] + 2. * k2['ode' ] + 2. * k3['ode' ] + k4['ode' ]) * (h / 6.)
    qf =      (k1['quad'] + 2. * k2['quad'] + 2. * k3['quad'] + k4['quad']) * (h / 6.)
    cf = 0.5 * (cs.mtimes(k1['r'].T, k1['r']) + 2. * cs.mtimes(k2['r'].T, k2['r']) + 2. * cs.mtimes(k3['r'].T, k3['r']) + cs.mtimes(k4['r'].T, k4['r'])) * (h / 6.)
    
    return cs.Function(name, [x0, u, t0], [xf, qf, cf], ['x0', 'p', 't0'], ['xf', 'qf', 'cf'])

t     = cs.MX.sym('t')       # time
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

#integrator = cs.integrator('pendulum_integrator', 'cvodes', {'x' : x, 'p' : u, 'ode' : dx, 'quad' : q}, {'tf' : ts})
integrator = rk4gn('pendulum_integrator', {'x' : x, 'p' : u, 'ode' : dx, 'quad' : q, 'r' : r}, {'tf' : ts})

integrator_out = integrator(x0 = x, p = u, t0 = t)
x_plus = integrator_out['xf']
qf     = integrator_out['qf']
cf     = integrator_out['cf']
integrator_sens = cs.Function('Integrator_Sensitivities', 
                              [t, x, u], 
                              [x_plus, cs.jacobian(x_plus, x), cs.jacobian(x_plus, u), 
                               qf    , cs.jacobian(    qf, x), cs.jacobian(    qf, u),
                               cf    , cs.jacobian(    cf, x), cs.jacobian(    cf, u)],
                              ['t0', 'x0', 'u'], 
                              ['xf', 'A', 'B', 'qf', 'qA', 'qB', 'cf', 'cA', 'cB'])

#------------------------------
# Generate C code for ODE model
#------------------------------
name = 'pendulum_ode'
#x_seed = cs.MX('x_seed', x.shape);
#u_seed = cs.MX('u_seed', u.shape);
ode = cs.Function(name, [t, x, u], 
                      [cs.densify(dx), cs.densify(cs.jacobian(dx, x)), cs.densify(cs.jacobian(dx, u)), 
                       cs.densify( q), cs.densify(cs.jacobian( q, x)), cs.densify(cs.jacobian( q, u)),
                       cs.densify( r), cs.densify(cs.jacobian( r, x)), cs.densify(cs.jacobian( r, u))], 
                      ['t', 'x0', 'u0'], ['xdot', 'A', 'B', 'q', 'qA', 'qB', 'r', 'rA', 'rB'])

x_seed = cs.MX.sym('x_seed', x.shape)
u_seed = cs.MX.sym('u_seed', u.shape)
ode_sens = cs.Function(name + "_sens", [t, x, u, x_seed, u_seed], 
                      [cs.densify(dx), cs.densify(cs.jtimes(dx, x, x_seed) + cs.jtimes(dx, u, u_seed))], 
                      ['t', 'x0', 'u0', 'x_seed', 'u_seed'], ['xdot', 'xdot_sens'])
gen = cs.CodeGenerator({'mex' : False, 'with_header' : True})
gen.add(ode)
gen.add(ode_sens)
name_c = '{0}_generated.c'.format(name)
name_h = '{0}_generated.h'.format(name)
gen.generate(name_c)
#os.rename(name_c, 'src/' + name_c)
#os.rename(name_h, 'src/' + name_h)

#------------------------------
# Generate test data
#------------------------------
x0 = cs.DM([1.0, 0.0])  # initial state

keys = ['t', 'x0', 'u', 'xdot', 'A_ode', 'B_ode', 'q', 'qA_ode', 'qB_ode', 
        'x_plus', 'A', 'B', 'qf', 'qA', 'qB', 'r', 'rA_ode', 'rB_ode', 
        'cf', 'cA', 'cB']
data = {}

for key in keys:
    data[key] = []

N = 600

for k in range(N):
    t_k = k * ts
    
    # input
    if t_k < 1.5:    
        u  = -10.0
    else:
        u = 0.0
                                        
    [x_plus, A, B, qf, qA, qB, cf, cA, cB] = integrator_sens(t_k, x0, u)
    [xdot, A_ode, B_ode, q, qA_ode, qB_ode, r, rA_ode, rB_ode] = ode(t_k, x0, u)
    
    data['t'     ].append(t_k)
    data['x0'    ].append(x0)
    data['u'     ].append(u)
    data['xdot'  ].append(xdot)
    data['A_ode' ].append(A_ode)
    data['B_ode' ].append(B_ode)
    data['q'     ].append(q)
    data['qA_ode'].append(qA_ode)
    data['qB_ode'].append(qB_ode)
    data['x_plus'].append(x_plus)
    data['A'     ].append(A)
    data['B'     ].append(B)
    data['qf'    ].append(qf)
    data['qA'    ].append(qA)
    data['qB'    ].append(qB)
    data['r'     ].append(r)
    data['rA_ode'].append(rA_ode)
    data['rB_ode'].append(rB_ode)
    data['cf'    ].append(cf)
    data['cA'    ].append(cA)
    data['cB'    ].append(cB)
    x0 = x_plus
    
ensure_dir_exist('data/rk4')
    
with open('data/rk4/pendulum.txt', 'w') as file:
    sep = ' '
    for k in range(N):
        for key in keys:
            np.array(data[key][k]).tofile(file, sep)        
            file.write('\n')
        
plt.subplot(2, 1, 1)
plt.step(cs.vertcat(data['t']), cs.transpose(cs.horzcat(*data['u' ]))      )
plt.subplot(2, 1, 2)
plt.plot(cs.vertcat(data['t']), cs.transpose(cs.horzcat(*data['x0'])), '.-')
plt.show()