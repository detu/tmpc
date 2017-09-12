#pragma once

#include "../qp/diagnostics.hpp"
#include "Trajectory.hpp"
#include "../qp/QpSize.hpp"
#include "../Matrix.hpp"
#include "../mpc/MpcQpSize.hpp"

#include <stdexcept>
#include <limits>
#include <iterator>
//#include <type_traits>

namespace tmpc
{
	/**
	 * \brief Vector of QpSize corresponding to an MPC problem with given sizes.
	 */
	std::vector<QpSize> RtiQpSize(std::size_t nt, std::size_t nx, std::size_t nu, std::size_t nc, std::size_t nct);

	/**
	 * \brief Implements an MPC controller with realtime iteration scheme.
	 *
	 * \tparam <D> Class defining problem dimensions
	 */
	template <typename Real_, typename D, typename OCP, typename QpWorkspace_>
	class RealtimeIteration
	{
		static auto constexpr NX = D::NX;
		static auto constexpr NU = D::NU;
		static auto constexpr NC = D::NC;
		static auto constexpr NCT = D::NCT;

	public:
		using Real = Real_;
		typedef QpWorkspace_ QpWorkspace;

		// TODO: parameterize Trajectory with K
		typedef Trajectory<NX, NU> WorkingPoint;

		RealtimeIteration(OCP const& ocp, WorkingPoint const& working_point)
		:	_ocp(ocp)
		,	workingPoint_(working_point)
		,	lowerBound_(working_point.nT(),
				constant<StaticVector<Real, NX>>(-inf_),
				constant<StaticVector<Real, NU>>(-inf_))
		,	upperBound_(working_point.nT(),
				constant<StaticVector<Real, NX>>( inf_),
				constant<StaticVector<Real, NU>>( inf_))
		,	qp_(mpcQpSize(working_point.nT(), NX, NU, NC, NCT))
		,	_prepared(false)
		{
			InitQP();
		}

		RealtimeIteration(RealtimeIteration const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		std::size_t nT() const { return workingPoint_.nT(); }

		// Feed current state x0, get back control input u.
		template <typename StateVector>
		decltype(auto) Feedback(const StateVector& x0)
		{
			if (!_prepared)
				throw std::logic_error("RealtimeIteration::Feedback(): RTI is not prepared.");

			/** embed current initial value */
			auto const w0 = eval(x0 - workingPoint_.get_x(0));
			qp_.problem().front().lbx(w0);
			qp_.problem().front().ubx(w0);

			/** solve QP */
			qp_.solve();

			// Add QP step to the working point.
			{
				for (std::size_t i = 0; i < nT() + 1; ++i)
					workingPoint_.set_x(i, workingPoint_.get_x(i) + qp_.solution()[i].x());

				for (std::size_t i = 0; i < nT(); ++i)
					workingPoint_.set_u(i, workingPoint_.get_u(i) + qp_.solution()[i].u());
			}

			_prepared = false;

			// Return the calculated control input.
			return workingPoint_.get_u(0);
		}

		/// \brief Number of iterations performed by the QP solver during the last Feedback phase.
		unsigned getNumIter() const { return qp_.numIter(); }

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("RealtimeIteration::Preparation(): RTI is already prepared.");

			// MK: Working point shifting is disabled. The main reason is that it is possible to have
			// time step in MPC different (typically bigger than) from the controller sampling time.
			// Perhaps, shifting policy can be defined externally.

			// shift(work_.workingPoint_);

			// Calculate new QP.
			UpdateQP();

			_prepared = true;
		}

		WorkingPoint const& getWorkingPoint() const { return workingPoint_; }

		void setWorkingPoint(WorkingPoint const& wp)
		{
			if (wp.nT() != workingPoint_.nT())
			{
				throw std::invalid_argument(
					"RealtimeIteration::setWorkingPoint(): new working point length is different. "
					"Changing prediction horizon on the fly is not supported (yet?).");
			}

			workingPoint_ = wp;
			_prepared = false;
		}

		// TODO: change LevenbergMarquardt to be an argument of Preparation()?
		Real getLevenbergMarquardt() const { return levenbergMarquardt_; }

		void setLevenbergMarquardt(Real val)
		{
			if (_prepared)
				throw std::logic_error("RealtimeIteration::setLevenbergMarquardt(): RTI must not be prepared when setting LevenbergMarquardt.");

			if (val < 0.)
				throw std::invalid_argument("RealtimeIteration::setLevenbergMarquardt(): LevenbergMarquardt must be non-negative.");

			if (levenbergMarquardt_ != val)
			{
				levenbergMarquardt_ = val;
				InitQP();	// Re-initialize QPs because of the changed levenbergMarquardt_ value.
			}
		}

		auto getQP() const
		{
			return qp_.problem();
		}

		auto getSolution() const
		{
			return qp_.solution();
		}

		// For setting options only!
		QpWorkspace& qpWorkspace()
		{
			return qp_;
		}

		/*
		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		unsigned nU() const	noexcept { return _Nu; }
		unsigned nX() const	noexcept { return _Nx; }
		unsigned nZ() const noexcept { return _Nz; }
		*/

	private:
		/**
		 * \brief Interface for setting MPC stage data.
		 * OCP::InitStage() and OCP::UpdateStage() use this interface to set Hessian of the stage cost,
		 * gradient of the stage cost, linearized stage constraints, state and input bounds.
		 *
		 * Hessian and gradient of Lagrange term:
		 * H = [Q, S
		 *      S, R]
		 * g = [q
		 *      r]
		 *
		 * NOTE:
		 * Potential issue with this implementation: it is easy to not to set one of the properties.
		 * But a single function with 11 arguments would look ugly.
		 * On the other hand, re-evaluation of the stage may be not necessary for some problems,
		 * and previous values would be kept. This allows more efficient implementations in some cases.
		 */
		class Stage
		{
		public:
			Stage(std::size_t i, RealtimeIteration& rti)
			:	i_(i)
			, 	work_(rti)
			,	qp_(rti.qp_.problem()[i])
			{
			}

			Stage(Stage const&) = delete;

			/**
			 * \brief On destruction, updates bounds of the QP.
			 */
			~Stage()
			{
				qp_.lbx(work_.lowerBound_.get_x(i_) - get_x());
				qp_.ubx(work_.upperBound_.get_x(i_) - get_x());
				qp_.lbu(work_.lowerBound_.get_u(i_) - get_u());
				qp_.ubu(work_.upperBound_.get_u(i_) - get_u());
			}

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return i_; }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return work_.workingPoint_.get_x(i_); }

			/**
			 * \brief Returns stage input.
			 */
			decltype(auto) get_u() const { return work_.workingPoint_.get_u(i_); }

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				qp_.Q(Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(Matrix const &R) {
				qp_.R(R + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(size_t i, size_t j, Matrix const &R)
			{
				auto constexpr M = rows<Matrix>();
				auto constexpr N = cols<Matrix>();

				auto const I = identity<NU, NU>();
				qp_.R(i, j, R + work_.levenbergMarquardt_ * submatrix(I, i, j, M, N));
										   // ^ Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_S(Matrix const &S) { qp_.S(S); }

			template <typename Matrix>
			void set_S(size_t i, size_t j, Matrix const &S)
			{
				qp_.S(i, j, S);
			}

			template <typename Vector>
			void set_q(Vector const &q) { qp_.q(q); }

			template <typename Vector>
			void set_r(Vector const &r) { qp_.r(r); }

			template <typename Vector>
			void set_r(size_t i, Vector const &r)
			{
				qp_.r(i, r);
			}

			template <typename Matrix>
			void set_A(Matrix const &A) { qp_.A(A); }

			template <typename Matrix>
			void set_B(Matrix const &B) { qp_.B(B); }

			template <typename Matrix>
			void set_B(size_t i, size_t j, Matrix const &B)
			{
				qp_.B(i, j, B);
			}

			template <typename Vector>
			void set_x_next(Vector const& x_next)
			{
				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				qp_.b(x_next - work_.workingPoint_.get_x(i_ + 1));
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { qp_.C(C); }

			template <typename Matrix>
			void set_D(Matrix const &val) { qp_.D(val); }

			template <typename Vector>
			void set_d_min(Vector const& d_min) { qp_.lbd(d_min); }

			template <typename Vector>
			void set_d_max(Vector const& d_max) { qp_.ubd(d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { work_.lowerBound_.set_x(i_, x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { work_.upperBound_.set_x(i_, x_max); }

			template <typename Vector>
			void set_u_min(Vector const& u_min)
			{
				work_.lowerBound_.set_u(i_, u_min);
			}

			template <typename Vector>
			void set_u_min(size_t i, Vector const& u_min)
			{
				work_.lowerBound_.set_u(i_, i, u_min);
			}

			template <typename Vector>
			void set_u_max(Vector const& u_max)
			{
				work_.upperBound_.set_u(i_, u_max);
			}

			template <typename Vector>
			void set_u_max(size_t i, Vector const& u_max)
			{
				work_.upperBound_.set_u(i_, i, u_max);
			}

		private:
			std::size_t const i_;
			RealtimeIteration& work_;
			typename QpWorkspace::Stage& qp_;
		};

		/**
		 * \brief Interface for setting MPC terminal stage data.
		 * OCP::InitTerminalStage() and OCP::UpdateTerminalStage() use this interface to set Hessian of the terminal stage cost,
		 * gradient of the terminal stage cost, linearized terminal stage constraints, and terminal state bounds.
		 *
		 * Hessian and gradient of Mayer term:
		 * H = [Q]
		 * g = [q]
		 *
		 */
		class TerminalStage
		{
		public:
			TerminalStage(RealtimeIteration& rti)
			:	work_(rti)
			,	qp_(rti.qp_.problem().back())
			{
			}

			TerminalStage(Stage const&) = delete;

			/**
			 * \brief On destruction, updates bounds of the QP.
			 */
			~TerminalStage()
			{
				qp_.lbx(work_.lowerBound_.get_x(get_index()) - get_x());
				qp_.ubx(work_.upperBound_.get_x(get_index()) - get_x());
			}

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return work_.workingPoint_.nT(); }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return work_.workingPoint_.get_x(get_index()); }

			// TODO: Isn't it more flexible and clear to add Lev-Mar in the OCP?
			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				qp_.Q(Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			decltype(auto) get_Q() const
			{
				return qp_.Q();
			}

			template <typename Vector>
			void set_q(Vector const &q)
			{
				qp_.q(q);
			}

			decltype(auto) get_q() const
			{
				return qp_.q();
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { qp_.C(C); }

			template <typename Vector>
			void set_d_end_min(Vector const& d_min) { qp_.lbd(d_min); }

			template <typename Vector>
			void set_d_end_max(Vector const& d_max) { qp_.ubd(d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { work_.lowerBound_.set_x(get_index(), x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { work_.upperBound_.set_x(get_index(), x_max); }

		private:
			RealtimeIteration& work_;
			typename QpWorkspace::Stage& qp_;
		};

		/**
		 * \brief Set all data of the QP to NaN.
		 */
		void SetQpNaN()
		{
			auto constexpr nan = std::numeric_limits<Real>::signaling_NaN();

			for (auto& qp : qp_.problem())
			{
				qp.Q(nan);
				qp.R(nan);
				qp.S(nan);
				qp.q(nan);
				qp.r(nan);

				qp.A(nan);
				qp.B(nan);
				qp.b(nan);

				qp.lbx(nan);
				qp.ubx(nan);
				qp.lbu(nan);
				qp.ubu(nan);

				qp.C(nan);
				qp.D(nan);

				qp.lbd(nan);
				qp.ubd(nan);
			}
		}

		// Uses ocp_ to initialize _QP based on current working point workingPoint_.
		void InitQP()
		{
			SetQpNaN();

			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, *this);
				_ocp.InitStage(stage);
			}

			TerminalStage terminal_stage(*this);
			_ocp.InitTerminalStage(terminal_stage);
		}

		// Uses ocp_ to update _QP based on current working point workingPoint_.
		void UpdateQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, *this);
				_ocp.UpdateStage(stage);
			}

			TerminalStage terminal_stage(*this);
			_ocp.UpdateTerminalStage(terminal_stage);

#if 0 && !defined(NDEBUG)
			// Error catching: check if we have any unexpected NaNs of infs in QP after update.

			std::vector<std::string> messages;
			qp::diagnose(qp_.problem(), std::back_inserter(messages));

			if (!messages.empty())
			{
				std::stringstream err_msg;
				err_msg << "RealtimeIteration encountered an ill-formed QP. Messages follow:" << std::endl;
				std::copy(messages.begin(), messages.end(), std::ostream_iterator<std::string>(err_msg, "\n"));

				throw std::logic_error(err_msg.str());
			}
#endif
		}

		// Private data members.
		OCP const& _ocp;	// <-- TODO: change to value semantics

		/// \brief Levenberg-Marquardt term.		
		Real levenbergMarquardt_ = 0.;

		/// \brief The working point (linearization point).
		WorkingPoint workingPoint_;

		/// \brief Lower bound of the working point.
		WorkingPoint lowerBound_;

		/// \brief Upper bound of the working point.
		WorkingPoint upperBound_;

		/// Multistage QP workspace.
		QpWorkspace qp_;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;

		static auto constexpr inf_ = std::numeric_limits<Real>::infinity();
	};
}
