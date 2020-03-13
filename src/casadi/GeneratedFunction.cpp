#include <tmpc/casadi/GeneratedFunction.hpp>


namespace tmpc :: casadi
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

        // Read In and Out argument sparsity
        size_t const n_in = fun_.f_->n_in();
        sparsityIn_.reserve(n_in);
        for (size_t i = 0; i < n_in; ++i)
            sparsityIn_.emplace_back(fun_.f_->sparsity_in(i));

        size_t const n_out = fun_.f_->n_out();
        sparsityOut_.reserve(n_out);
        for (size_t i = 0; i < n_out; ++i)
            sparsityOut_.emplace_back(fun_.f_->sparsity_out(i));

        allocateDataInOut();
    }


    GeneratedFunction::GeneratedFunction(GeneratedFunction const& rhs)
    :	fun_ {rhs.fun_}
    ,	name_ {rhs.name_}
    ,   sparsityIn_ {rhs.sparsityIn_}
    ,   sparsityOut_ {rhs.sparsityOut_}
    ,	arg_(rhs.arg_.size())
    ,	res_(rhs.res_.size())
    ,	iw_(rhs.iw_.size())
    ,	w_(rhs.w_.size())
    {
        allocateDataInOut();
    }


    GeneratedFunction::~GeneratedFunction()
    {
    }


    void GeneratedFunction::allocateDataInOut()
    {
        auto const n_in = this->n_in();
        auto const n_out = this->n_out();

        dataIn_.reserve(n_in);
        for (size_t i = 0; i < n_in; ++i)
            dataIn_.emplace_back(new casadi_real[sparsityIn_[i].nnz()]);

        dataOut_.reserve(n_out);
        for (size_t i = 0; i < n_out; ++i)
            dataOut_.emplace_back(new casadi_real[sparsityOut_[i].nnz()]);
    }
}