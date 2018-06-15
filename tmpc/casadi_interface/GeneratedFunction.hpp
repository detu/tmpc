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
#include <sstream>

#include <casadi/mem.h>


namespace casadi_interface
{
	/// TODO: implement move constructor.
	class GeneratedFunction
	{
	public:
		/// @brief Constructor.
		GeneratedFunction(casadi_functions const * f, std::string const& n = "Unknown CasADi function")
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


		/// @brief Copy constructor.
		GeneratedFunction(GeneratedFunction const& rhs)
		:	fun_{rhs.fun_}
		,	name_{rhs.name_}
		,	arg_(rhs.arg_.size())
		,	res_(rhs.res_.size())
		,	iw_(rhs.iw_.size())
		,	w_(rhs.w_.size())
		{
		}
		

		// MK: I am too lazy to implement the assignment operator at the moment, 
		// so I prevent it from being used.
		// TODO: implement assignment and move-assignment.
		GeneratedFunction& operator=(GeneratedFunction const& rhs) = delete;


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

		void operator()(std::initializer_list<const casadi_real *> arg, std::initializer_list<casadi_real *> res) const
		{
			if (arg.size() != n_in())
				throw std::invalid_argument("Invalid number of input arguments to " + name());

			if (res.size() != n_out())
				throw std::invalid_argument("Invalid number of output arguments to " + name());

			std::copy(arg.begin(), arg.end(), arg_.begin());
			std::copy(res.begin(), res.end(), res_.begin());

			fun_.f_->eval(arg_.data(), res_.data(), iw_.data(), w_.data(), 0);
		}

	private:
		struct Functions
		{
			Functions(casadi_functions const * f) 
			:	f_(f)
			{ 
				f_->incref(); 
			}


			Functions(Functions const& f) 
			:	f_(f.f_)
			{ 
				f_->incref();
			}


			Functions& operator=(Functions const& rhs) noexcept
			{
				if (&rhs != this)
				{
					f_->decref();
					f_ = rhs.f_;
					f_->incref();
				}

				return *this;
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
		mutable std::vector<casadi_real const *> arg_;
		mutable std::vector<casadi_real *> res_;
		mutable std::vector<casadi_int> iw_;
		mutable std::vector<casadi_real> w_;
	};

	// Some interesting ideas here: https://habrahabr.ru/post/228031/
	/*
	template <typename... Args>
	void call(GeneratedFunction const& f, Args&&... args)
	{
		std::tuple<Args...> args_ = std::make_tuple(std::forward<Args>(args)...);
	}
	*/
}
