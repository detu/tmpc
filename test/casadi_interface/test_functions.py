import casadi as cs
        
A = cs.MX.sym('A', 3, 2)
B = cs.MX.sym('B', 2, 2)
x = cs.MX.sym('x'      )

f = cs.Function('f', [A, B, x], [cs.mtimes(A * x, B), cs.sum1(cs.mtimes(A, B))])

def functions():
    return [f]
