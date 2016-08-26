clear all;clc;
Ts = 0.1;

% Differential states, controls and online data
DifferentialState x1 x2
Control u
OnlineData thetaR1 kappaR1 thetaR11
x = diffStates;

NX = length(diffStates);
NU = length(controls);

% Chen1998 dynamics
mu = 0.5;
f = dot(x) == [x2 + u*(mu+(1-mu)*x1); x1 + u*(mu-4*(1-mu)*x2)];

% Stage cost (least-squares vector)
h = [x1;x2;u];

% Terminal cost (least-squares vector)
hN = [x1;x2];

% Definition of external cost matrices
NY = length(h);
NYN = length(hN);
W = diag([0.5,0.5,1.0]);
WN = [16.5926, 11.5926; 11.5926, 16.5926];
% W = acado.DMatrix(W_mat);
% WN = acado.DMatrix(WN_mat);

N = 13;
ocp = acado.OCP( 0.0, N*Ts, N );

% Define OCP objective
ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );

% Inequality constraints
ocp.subjectTo( -2 <= u <= 2 );
ocp.setModel(f);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',         'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',           'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',            'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',               'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',          N                   );
mpc.set( 'QP_SOLVER',                     'QP_QPOASES'        );
mpc.set( 'GENERATE_MAKE_FILE',            'NO'                );
mpc.set( 'GENERATE_TEST_FILE',            'NO'                );
mpc.set( 'GENERATE_MATLAB_INTERFACE',     'NO'                );
mpc.set( 'CG_USE_OPENMP',                 'NO'                );

mpc.exportCode( 'export_MPC' );
