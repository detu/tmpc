/*
 * CasADiGeneratedFunction.h
 *
 *  Created on: Apr 19, 2016
 *      Author: kotlyar
 */

#pragma once

#include <vector>
#include <string>
#include <array>
//#include <initializer_list>
#include <sstream>

#ifdef real_t
#undef real_t
#endif

#define CASADI_GENERATED_FUNCTION_INTERFACE(name, N_IN, N_OUT) /*#name,*/ name, name##_incref, name##_decref, name##_n_in, name##_n_out, \
	name##_name_in, name##_name_out, name##_sparsity_in, name##_sparsity_out, name##_work, N_IN, N_OUT
#define CASADI_GENERATED_FUNCTION_CLASS(name, N_IN, N_OUT) \
	::casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(name, N_IN, N_OUT)>

namespace casadi_interface
{
	typedef double real_t;

	template <
		//char const * Name,
		int (*_fun)(const real_t** arg, real_t** res, int* iw, real_t* w, int mem),
		void (*_fun_incref)(void),
		void (*_fun_decref)(void),
		int (*_fun_n_in)(void),
		int (*_fun_n_out)(void),
		const char* (*_fun_name_in)(int i),
		const char* (*_fun_name_out)(int i),
		const int* (*_fun_sparsity_in)(int i),
		const int* (*_fun_sparsity_out)(int i),
		int (*_fun_work)(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w),
		std::size_t N_IN,
		std::size_t N_OUT
	>
	class GeneratedFunction
	{
	public:
		GeneratedFunction()
		{
			// Check number of inputs and outputs.
			if (_fun_n_in() != N_IN)
				throw std::logic_error("Generated function's number of inputs is different from expected!");

			if (_fun_n_out() != N_OUT)
				throw std::logic_error("Generated function's number of outputs is different from expected!");

			// Get work memory size.
			int sz_arg, sz_res, sz_iw, sz_w;
			if (int code = _fun_work(&sz_arg, &sz_res, &sz_iw, &sz_w) != 0)
				throw std::runtime_error(name() + "_fun_work() returned " + std::to_string(code));

			_arg.resize(sz_arg);
			_res.resize(sz_res);
			_iw.resize(sz_iw);
			_w.resize(sz_w);
		}

		static std::string const& name()
		{
			static std::string const n = "[FUNCTION NAME]";
			return n;
		}

		static std::size_t constexpr n_in () { return N_IN ; }
		static std::size_t constexpr n_out() { return N_OUT; }

		int n_row_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				throw std::out_of_range("GeneratedFunction::n_row_in(): index is out of range");

			return _fun_sparsity_in(ind)[0];
		}

		int n_col_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				throw std::out_of_range("GeneratedFunction::n_col_in(): index is out of range");

			return _fun_sparsity_in(ind)[1];
		}

		int n_row_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				throw std::out_of_range("GeneratedFunction::n_row_out(): index is out of range");

			return _fun_sparsity_out(ind)[0];
		}

		int n_col_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				throw std::out_of_range("GeneratedFunction::n_col_out(): index is out of range");

			return _fun_sparsity_out(ind)[1];
		}

		void operator()(std::array<const real_t *, N_IN> arg, std::array<real_t *, N_OUT> res) const
		{
			std::copy(arg.begin(), arg.end(), _arg.begin());
			std::copy(res.begin(), res.end(), _res.begin());

			_fun(_arg.data(), _res.data(), _iw.data(), _w.data(), 0);
		}

	private:
		struct RefHolder
		{
			RefHolder() { _fun_incref(); }
			~RefHolder() { _fun_decref(); }
		};

		RefHolder refHolder_;

		//static std::string const _name { "Name" };

		/*
		"To allow the evaluation to be performed efficiently with a small memory footprint, the
		user is expected to pass four work arrays. The function fname_work returns the length
		of these arrays, which have entries of type const double*, double*, int and double,
		respectively."

		CasADi user guide, section 5.3.
		http://casadi.sourceforge.net/users_guide/casadi-users_guide.pdf
		*/
		mutable std::vector<real_t const *> _arg;
		mutable std::vector<real_t *> _res;
		mutable std::vector<int> _iw;
		mutable std::vector<real_t> _w;
	};

	// Some interesting ideas here: https://habrahabr.ru/post/228031/
	/*
	template <typename... Args>
	void call(GeneratedFunction const& f, Args&&... args)
	{
		std::tuple<Args...> args_ = std::make_tuple(std::forward<Args>(args)...);
	}
	*/
} /* namespace mpmc */
