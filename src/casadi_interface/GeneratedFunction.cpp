#include <tmpc/casadi_interface/GeneratedFunction.hpp>


namespace casadi_interface
{
    GeneratedFunction::GeneratedFunction(casadi_functions const * f, std::string const& n)
    :	fun_(f)
    ,	name_(n)
    {
        // Get work memory size.
        casadi_int sz_arg, sz_res, sz_iw, sz_w;
        if (int code = fun_.f_->work(&sz_arg, &sz_res, &sz_iw, &sz_w) != 0)
            throw std::runtime_error(name() + "_fun_work() returned " + std::to_string(code));

        arg_.resize(sz_arg);
        res_.resize(sz_res);
        iw_.resize(sz_iw);
        w_.resize(sz_w);
    }


    GeneratedFunction::GeneratedFunction(GeneratedFunction const& rhs)
    :	fun_{rhs.fun_}
    ,	name_{rhs.name_}
    ,	arg_(rhs.arg_.size())
    ,	res_(rhs.res_.size())
    ,	iw_(rhs.iw_.size())
    ,	w_(rhs.w_.size())
    {
    }


    GeneratedFunction::~GeneratedFunction()
    {        
    }
}