import casadi as cs
import sys
import os.path
import importlib.util


#------------------------------
# Import module
#------------------------------
module_path = sys.argv[1]
module_file_name = os.path.basename(module_path)
name_c = os.path.splitext(module_file_name)[0] + ".c"

print("Importing {0} ({1})...".format(name_c, module_path))
spec = importlib.util.spec_from_file_location("module", module_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

#------------------------------
# Generate C code
#------------------------------
funcs = module.functions()
print("Found {0} functions".format(len(funcs)))

gen = cs.CodeGenerator(name_c, {'mex' : False, 'with_header' : True, 'with_mem' : True})
for f in funcs:
    gen.add(f)

gen.generate()
#input('PRESS ANY KEY...')