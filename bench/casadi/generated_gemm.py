import casadi as cs


def generated_gemm(m, n, k, type_name):
    if type_name == 'MX':
        T = cs.MX
    elif type_name == 'SX':
        T = cs.SX
    else:
        raise ValueError('Invalid value of type_name')

    A = T.sym('A', k, m)
    B = T.sym('B', k, n)
    C = T.sym('C', m, n)

    return cs.Function('generated_gemm_{0}x{1}x{2}_{3}'.format(m, n, k, type_name), [A, B, C], [C + cs.mtimes(cs.transpose(A), B)])


def functions():
    f = []
    for T in ['MX', 'SX']:
        for m in range(1, 41):
            f.append(generated_gemm(m, m, m, T))
    
    return f