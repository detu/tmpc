/*
 * CasADiGeneratedFunction.h
 *
 *  Created on: Apr 19, 2016
 *      Author: kotlyar
 */

#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/casadi/Sparsity.hpp>
#include <tmpc/casadi/CompressedColumnStorage.hpp>
#include <tmpc/Exception.hpp>

#include <vector>
#include <string>
#include <array>
#include <initializer_list>
#include <sstream>
#include <tuple>
#include <utility>
#include <memory>

#include <casadi/mem.h>

#include <blaze/Math.h>


namespace tmpc :: casadi
{
	/// TODO: implement move constructor.
	class GeneratedFunction
	{
	public:
		/// @brief Constructor.
		GeneratedFunction(casadi_functions const * f, std::string const& n = "CasADi function");

		/// @brief Copy constructor.
		GeneratedFunction(GeneratedFunction const& rhs);

		~GeneratedFunction();
		

		// MK: I am too lazy to implement the assignment operator at the moment, 
		// so I prevent it from being used.
		// TODO: implement assignment and move-assignment.
		GeneratedFunction& operator=(GeneratedFunction const& rhs) = delete;


		std::string const& name() const noexcept
		{
			return name_;
		}


		size_t const n_in () const noexcept { return size(sparsityIn_); }
		size_t const n_out() const noexcept { return size(sparsityOut_); }

		
		int n_row_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				TMPC_THROW_EXCEPTION(std::out_of_range("input argument index out of range"));

			return sparsityIn_[ind].rows();
		}


		int n_col_in(int ind) const
		{
			if (ind < 0 || ind >= n_in())
				TMPC_THROW_EXCEPTION(std::out_of_range("input argument index out of range"));

			return sparsityIn_[ind].columns();
		}


		int n_row_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				TMPC_THROW_EXCEPTION(std::out_of_range("output argument index out of range"));

			return sparsityOut_[ind].rows();
		}


		int n_col_out(int ind) const
		{
			if (ind < 0 || ind >= n_out())
				TMPC_THROW_EXCEPTION(std::out_of_range("output argument index out of range"));

			return sparsityOut_[ind].columns();
		}


		void operator()(std::initializer_list<const casadi_real *> arg, std::initializer_list<casadi_real *> res) const
		{
			if (arg.size() != n_in())
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid number of input arguments to " + name()));

			if (res.size() != n_out())
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid number of output arguments to " + name()));

			std::copy(arg.begin(), arg.end(), arg_.begin());
			std::copy(res.begin(), res.end(), res_.begin());

			fun_.f_->eval(arg_.data(), res_.data(), iw_.data(), w_.data(), 0);
		}


		template <typename... ArgsIn, typename... ArgsOut>
		void operator()(std::tuple<ArgsIn&...> in, std::tuple<ArgsOut&...> out) const
		{
			// Check number of input arguments
			if (sizeof...(ArgsIn) != n_in())
				TMPC_THROW_EXCEPTION(std::invalid_argument("Wrong number of input arguments"));

			// Check number of output arguments
			if (sizeof...(ArgsOut) != n_out())
				TMPC_THROW_EXCEPTION(std::invalid_argument("Wrong number of output arguments"));

			// Copy data from input arguments to dataIn_ arrays 
			// and/or set pointers, depending on the argument types.
			copyIn(std::index_sequence_for<ArgsIn...>(), in);

			// Set pointers to output arguments, depending on the argument types.
			setOutPtr(std::index_sequence_for<ArgsOut...>(), out);

			// Evaluate the function
			fun_.f_->eval(arg_.data(), res_.data(), iw_.data(), w_.data(), 0);

			// Copy data from dataOut_ arrays to output arguments if neccessary,
			// depending on argument types.
			copyOut(std::index_sequence_for<ArgsOut...>(), out);
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
		std::vector<Sparsity> sparsityIn_;
		std::vector<Sparsity> sparsityOut_;

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
		mutable std::vector<std::unique_ptr<casadi_real[]>> dataIn_;
		mutable std::vector<std::unique_ptr<casadi_real[]>> dataOut_;


		template <std::size_t... Is, typename... Args>
		void copyIn(std::index_sequence<Is...>, std::tuple<Args&...> args) const
		{
			(copyIn(Is, std::get<Is>(args)), ...);
		}


		template <typename MT, bool SO>
		void copyIn(size_t i, blaze::Matrix<MT, SO> const& arg) const
		{
			toCompressedColumnStorage(arg, dataIn_[i].get(), sparsityIn_[i]);
			arg_[i] = dataIn_[i].get();
		}


		template <typename VT, bool TF>
		void copyIn(size_t i, blaze::Vector<VT, TF> const& arg) const
		{
			toCompressedColumnStorage(arg, dataIn_[i].get(), sparsityIn_[i]);
			arg_[i] = dataIn_[i].get();
		}


		void copyIn(size_t i, casadi_real const& arg) const
		{
			if (sparsityIn_[i].nnz() != 1)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid size of input argument " + std::to_string(i) + " in CasADi function " + name()));

			arg_[i] = &arg;
		}


		void copyIn(size_t i, std::nullptr_t) const
		{
			arg_[i] = nullptr;
		}


		template <std::size_t... Is, typename... Args>
		void setOutPtr(std::index_sequence<Is...>, std::tuple<Args&...> args) const
		{
			(setOutPtr(Is, std::get<Is>(args)), ...);
		}


		template <typename MT, bool SO>
		void setOutPtr(size_t i, blaze::Matrix<MT, SO>& arg) const
		{
			res_[i] = dataOut_[i].get();
		}


		template <typename VT, bool TF>
		void setOutPtr(size_t i, blaze::Vector<VT, TF>& arg) const
		{
			res_[i] = dataOut_[i].get();
		}


		void setOutPtr(size_t i, casadi_real& arg) const
		{
			res_[i] = &arg;
		}


		void setOutPtr(size_t i, std::nullptr_t) const
		{
			res_[i] = nullptr;
		}


		template <std::size_t... Is, typename... Args>
		void copyOut(std::index_sequence<Is...>, std::tuple<Args&...> args) const
		{
			(copyOut(Is, std::get<Is>(args)), ...);
		}


		template <typename MT, bool SO>
		void copyOut(size_t i, blaze::Matrix<MT, SO>& arg) const
		{
			fromCompressedColumnStorage(res_[i], sparsityOut_[i], arg);
		}


		template <typename VT, bool TF>
		void copyOut(size_t i, blaze::Vector<VT, TF>& arg) const
		{
			fromCompressedColumnStorage(res_[i], sparsityOut_[i], arg);
		}


		void copyOut(size_t i, casadi_real& arg) const
		{
			arg = *res_[i];
		}


		void copyOut(size_t i, std::nullptr_t) const
		{
		}


		void allocateDataInOut();
	};
}
