import casadi as cs
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
        
A = cs.MX.sym('A', 3, 2)
B = cs.MX.sym('B', 2, 2)
x = cs.MX.sym('x'      )

f = cs.Function('f', [A, B, x], [cs.mtimes(A * x, B), cs.sum1(cs.mtimes(A, B))])

#------------------------------
# Generate C code
#------------------------------
name_c = sys.argv[1]
gen = cs.CodeGenerator(name_c, {'mex' : False, 'with_header' : True})
gen.add(f)
gen.generate()
