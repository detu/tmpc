/*
 * CasADiGeneratedFunction.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: kotlyar
 */

#include <CasADiGeneratedFunction.hpp>

#include <stdexcept>
#include <sstream>

namespace mpmc
{
	CasADiGeneratedFunction::CasADiGeneratedFunction(
		const std::string& name,
		int (*fun)(const real_t** arg, real_t** res, int* iw, real_t* w, int mem),
		void (*fun_incref)(void),
		void (*fun_decref)(void),
		int (*fun_n_in)(void),
		int (*fun_n_out)(void),
		const char* (*fun_name_in)(int i),
		const char* (*fun_name_out)(int i),
		const int* (*fun_sparsity_in)(int i),
		const int* (*fun_sparsity_out)(int i),
		int (*fun_work)(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w))
	:	_name(name),
		_fun(fun),
		_fun_incref(fun_incref),
		_fun_decref(fun_decref),
		_fun_n_in(fun_n_in),
		_fun_n_out(fun_n_out),
		_fun_name_in(fun_name_in),
		_fun_name_out(fun_name_out),
		_fun_sparsity_in(fun_sparsity_in),
		_fun_sparsity_out(fun_sparsity_out),
		_fun_work(fun_work)
	{
		_fun_incref();

		int sz_arg, sz_res, sz_iw, sz_w;
		if (int code = _fun_work(&sz_arg, &sz_res, &sz_iw, &sz_w) != 0)
			throw std::runtime_error(name + "_fun_work() returned " + std::to_string(code));

		_arg.resize(sz_arg);
		_res.resize(sz_res);
		_iw.resize(sz_iw);
		_w.resize(sz_w);
	}

	CasADiGeneratedFunction::~CasADiGeneratedFunction()
	{
		_fun_decref();
	}

	const std::string& CasADiGeneratedFunction::name() const
	{
		return _name;
	}

	int CasADiGeneratedFunction::n_in() const
	{
		return _fun_n_in();
	}

	int CasADiGeneratedFunction::n_out() const
	{
		return _fun_n_out();
	}

	void CasADiGeneratedFunction::operator()(std::initializer_list<const real_t *> arg, std::initializer_list<real_t *> res)
	{
		if (arg.size() != n_in())
		{
			std::stringstream msg;
			msg << "Invalid number of input arguments passed to CasADi function \"" << name() << "\". "
					<< "Expected " << n_in() << ", got " << arg.size() << ".";
			throw std::logic_error(msg.str());
		}

		if (res.size() != n_out())
		{
			std::stringstream msg;
			msg << "Invalid number of output arguments passed to CasADi function \"" << name() << "\". "
					<< "Expected " << n_out() << ", got " << res.size() << ".";
			throw std::logic_error(msg.str());
		}

		std::copy(arg.begin(), arg.end(), _arg.begin());
		std::copy(res.begin(), res.end(), _res.begin());

		_fun(_arg.data(), _res.data(), _iw.data(), _w.data(), 0);
	}
} /* namespace mpmc */
