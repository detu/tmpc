import casadi as cs


def _matrixMultiply(m, n, p, type_name):
    if type_name == 'MX':
        T = cs.MX
    elif type_name == 'SX':
        T = cs.SX
    else:
        raise ValueError('Invalid value of type_name')

    A = T.sym('A', m, n)
    B = T.sym('B', n, p)

    return cs.Function('MTimes_{0}x{1}x{2}_{3}'.format(m, n, p, type_name), [A, B], [cs.mtimes(A, B)])


def functions():
    f = []
    for T in ['MX', 'SX']:
        for m in [2, 3, 5, 10, 20, 30]:
            f.append(_matrixMultiply(m, m, m, T))
    
    return f