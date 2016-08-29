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
#include <initializer_list>

#ifdef real_t
#undef real_t
#endif

#define CASADI_GENERATED_FUNCTION_INTERFACE(name) #name, name, name##_incref, name##_decref, name##_n_in, name##_n_out, \
	name##_name_in, name##_name_out, name##_sparsity_in, name##_sparsity_out, name##_work

namespace casadi_interface
{
	class GeneratedFunction
	{
	public:
		typedef double real_t;

		GeneratedFunction(
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
			int (*fun_work)(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w));

		~GeneratedFunction();

		const std::string& name() const;
		int n_in() const;
		int n_out() const;
		int n_row_in(int ind) const;
		int n_col_in(int ind) const;
		int n_row_out(int ind) const;
		int n_col_out(int ind) const;

		//void operator()(std::array<const real_t *, 3> arg, std::array<real_t *, 3> res);
		void operator() (std::initializer_list<const real_t *> arg, std::initializer_list<real_t *> res) const;

	private:
		const std::string _name;
		int (*_fun)(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
		void (*_fun_incref)(void);
		void (*_fun_decref)(void);
		int (*_fun_n_in)(void);
		int (*_fun_n_out)(void);
		const char* (*_fun_name_in)(int i);
		const char* (*_fun_name_out)(int i);
		const int* (*_fun_sparsity_in)(int i);
		const int* (*_fun_sparsity_out)(int i);
		int (*_fun_work)(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);

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
