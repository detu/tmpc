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

#include <casadi/mem.h>

namespace casadi_interface
{
	class GeneratedFunction
	{
	public:
		GeneratedFunction(casadi_functions const * f, std::string const& n = "Unknown CasADi function")
		:	fun_(f)
		,	name_(n)
		{
			// Get work memory size.
			int sz_arg, sz_res, sz_iw, sz_w;
			if (int code = fun_.f_->work(&sz_arg, &sz_res, &sz_iw, &sz_w) != 0)
				throw std::runtime_error(name() + "_fun_work() returned " + std::to_string(code));

			_iw.resize(sz_iw);
			_w.resize(sz_w);
		}

		std::string const& name() const
		{
			return name_;
		}

		std::size_t const n_in () const { return fun_.f_->n_in() ; }
		std::size_t const n_out() const { return fun_.f_->n_out(); }

		int n_row_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				throw std::out_of_range("GeneratedFunction::n_row_in(): index is out of range");

			return fun_.f_->sparsity_in(ind)[0];
		}

		int n_col_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				throw std::out_of_range("GeneratedFunction::n_col_in(): index is out of range");

			return fun_.f_->sparsity_in(ind)[1];
		}

		int n_row_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				throw std::out_of_range("GeneratedFunction::n_row_out(): index is out of range");

			return fun_.f_->sparsity_out(ind)[0];
		}

		int n_col_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				throw std::out_of_range("GeneratedFunction::n_col_out(): index is out of range");

			return fun_.f_->sparsity_out(ind)[1];
		}

		void operator()(std::initializer_list<const casadi_real_t *> arg, std::initializer_list<casadi_real_t *> res) const
		{
			if (arg.size() != n_in())
				throw std::invalid_argument("Invalid number of input arguments to " + name());

			if (res.size() != n_out())
				throw std::invalid_argument("Invalid number of output arguments to " + name());

			// See https://github.com/casadi/casadi/issues/2091 for the explanation why the const_cast<>'s are here.
			fun_.f_->eval(const_cast<casadi_real_t const **>(arg.begin()), const_cast<casadi_real_t **>(res.begin()), _iw.data(), _w.data(), 0);
		}

	private:
		struct Functions
		{
			Functions(casadi_functions const * f) 
			:	f_(f)
			{ 
				f_->incref(); 
			}

			~Functions() 
			{ 
				f_->decref(); 
			}

			casadi_functions const * f_;
		};

		Functions const fun_;
		std::string const name_;

		//static std::string const _name { "Name" };

		/*
		"To allow the evaluation to be performed efficiently with a small memory footprint, the
		user is expected to pass four work arrays. The function fname_work returns the length
		of these arrays, which have entries of type const double*, double*, int and double,
		respectively."

		CasADi user guide, section 5.3.
		http://casadi.sourceforge.net/users_guide/casadi-users_guide.pdf
		*/
		mutable std::vector<int> _iw;
		mutable std::vector<casadi_real_t> _w;
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
