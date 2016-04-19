/*
 * CasADiGeneratedFunction.h
 *
 *  Created on: Apr 19, 2016
 *      Author: kotlyar
 */

#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <array>
#include <initializer_list>

#define CASADI_GENERATED_FUNCTION_INTERFACE(name) #name, name, name##_incref, name##_decref, name##_n_in, name##_n_out, \
	name##_name_in, name##_name_out, name##_sparsity_in, name##_sparsity_out, name##_work

namespace mpmc
{
	class CasADiGeneratedFunction
	{
	public:
		typedef double real_t;

		CasADiGeneratedFunction(
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
			int (*fun_work)(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w)) throw (std::runtime_error);

		~CasADiGeneratedFunction();

		const std::string& name() const;
		int n_in() const;
		int n_out() const;

		void operator()(std::array<const real_t *, 3> arg, std::array<real_t *, 3> res);
		//void operator()(std::initializer_list<const real_t *> arg, std::initializer_list<real_t *> res);

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
		std::vector<const real_t *> _arg;
		std::vector<real_t *> _res;
		std::vector<int> _iw;
		std::vector<real_t> _w;
	};
} /* namespace mpmc */
