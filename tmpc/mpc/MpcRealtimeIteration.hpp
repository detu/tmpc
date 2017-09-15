#pragma once

#include <tmpc/qp/diagnostics.hpp>
#include <tmpc/qp/QpStageBase.hpp>
#include <tmpc/mpc/MpcTrajectoryPoint.hpp>
#include <tmpc/mpc/MpcSize.hpp>

#include <stdexcept>
#include <limits>
#include <iterator>
#include <vector>

namespace tmpc
{
	/**
	 * \brief Implements an MPC controller with realtime iteration scheme.
	 *
	 * \tparam <Real_> Real type
	 */
	template <typename Real_, typename OCP, typename QpWorkspace_>
	class MpcRealtimeIteration
	{
	public:
		using Real = Real_;
		using QpWorkspace = QpWorkspace_;
		using TrajectoryPoint = MpcTrajectoryPoint<Real>;
		using WorkingPoint = std::vector<TrajectoryPoint>;

		MpcRealtimeIteration(OCP const& ocp, WorkingPoint const& working_point)
		:	ocp_(ocp)
		,	qp_(ocp.dimensions())
		,	workingPoint_(working_point)
		,	prepared_(false)
		{
			// Initialize bounds to +-inf
			lowerBound_.reserve(working_point.size());
			upperBound_.reserve(working_point.size());

			for (auto const& sz : ocp.dimensions())
			{
				lowerBound_.emplace_back(
					DynamicVector<Real>(sz.nx(), -inf()),
					DynamicVector<Real>(sz.nu(), -inf())
				);

				upperBound_.emplace_back(
					DynamicVector<Real>(sz.nx(), inf()),
					DynamicVector<Real>(sz.nu(), inf())
				);
			}

			// Initialize the QP at working point
			InitQP();
		}

		MpcRealtimeIteration(MpcRealtimeIteration const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		size_t nT() const { return workingPoint_.size() - 1; }

		// Feed current state x0, get back control input u.
		template <typename StateVector>
		decltype(auto) Feedback(const StateVector& x0)
		{
			if (!prepared_)
				throw std::logic_error("MpcRealtimeIteration::Feedback(): RTI is not prepared.");

			/** embed current initial value */
			auto const w0 = eval(x0 - workingPoint_[0].x());
			qp_.problem()[0].lbx(w0);
			qp_.problem()[0].ubx(w0);

			/** solve QP */
			qp_.solve();

			// Add the calculated QP step to the working point.
			for (size_t i = 0; i < nT() + 1; ++i)
			{
				workingPoint_[i].x(workingPoint_[i].x() + qp_.solution()[i].x());
				workingPoint_[i].u(workingPoint_[i].u() + qp_.solution()[i].u());
			}
			
			prepared_ = false;

			// Return the calculated control input.
			return workingPoint_[0].u();
		}

		/// \brief Number of iterations performed by the QP solver during the last Feedback phase.
		unsigned numIter() const { return qp_.numIter(); }

		void Preparation()
		{
			if (prepared_)
				throw std::logic_error("MpcRealtimeIteration::Preparation(): RTI is already prepared.");

			// MK: Working point shifting is disabled. The main reason is that it is possible to have
			// time step in MPC different (typically bigger than) from the controller sampling time.
			// Perhaps, shifting policy can be defined externally.

			// shift(rti_.workingPoint_);

			// Calculate new QP.
			UpdateQP();

			prepared_ = true;
		}

		WorkingPoint const& workingPoint() const { return workingPoint_; }

		void workingPoint(WorkingPoint const& wp)
		{
			if (wp.size() != workingPoint_.size())
			{
				throw std::invalid_argument(
					"MpcRealtimeIteration::workingPoint(): new working point length is different. "
					"Changing prediction horizon on the fly is not supported (yet?).");
			}

			workingPoint_ = wp;
			prepared_ = false;
		}

		// TODO: change LevenbergMarquardt to be an argument of Preparation()?
		Real levenbergMarquardt() const { return levenbergMarquardt_; }

		void levenbergMarquardt(Real val)
		{
			if (prepared_)
				throw std::logic_error("MpcRealtimeIteration::levenbergMarquardt(): RTI must not be prepared when setting LevenbergMarquardt.");

			if (val < 0.)
				throw std::invalid_argument("MpcRealtimeIteration::levenbergMarquardt(): LevenbergMarquardt must be non-negative.");

			if (levenbergMarquardt_ != val)
			{
				levenbergMarquardt_ = val;
				InitQP();	// Re-initialize QPs because of the changed levenbergMarquardt_ value.
			}
		}

		QpWorkspace const& qp() const
		{
			return qp_;
		}
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
			Stage(std::size_t i, MpcRealtimeIteration& rti)
			:	i_(i)
			,	rti_(rti)
			,	qp_(rti.qp_.problem()[i])
			,	workingPoint_(rti_.workingPoint_[i])
			,	lowerBound_(rti_.lowerBound_[i])
			,	upperBound_(rti_.upperBound_[i])
			{
			}

			Stage(Stage const&) = delete;

			/**
			 * \brief On destruction, updates bounds of the QP.
			 */
			~Stage()
			{
				qp_.lbx(lowerBound_.x() - workingPoint_.x());
				qp_.ubx(upperBound_.x() - workingPoint_.x());
				qp_.lbu(lowerBound_.u() - workingPoint_.u());
				qp_.ubu(upperBound_.u() - workingPoint_.u());
			}

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return i_; }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return workingPoint_.x(); }

			/**
			 * \brief Returns stage input.
			 */
			decltype(auto) get_u() const { return workingPoint_.u(); }

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				qp_.Q(Q + IdentityMatrix<Real>(qp_.size().nx()) * rti_.levenbergMarquardt_);	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(Matrix const &R) {
				qp_.R(R + IdentityMatrix<Real>(qp_.size().nu()) * rti_.levenbergMarquardt_);	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(size_t i, size_t j, Matrix const &R)
			{
				auto const M = rows(R);
				auto const N = columns(R);

				auto const I = IdentityMatrix<Real>(qp_.size().nu());
				qp_.R(i, j, R + rti_.levenbergMarquardt_ * submatrix(I, i, j, M, N));
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
				qp_.b(x_next - rti_.workingPoint_[i_ + 1].x());
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { qp_.C(C); }

			template <typename Matrix>
			void set_D(Matrix const &val) { qp_.D(val); }

			template <typename Vector>
			void set_d_min(Vector const& d_min) { qp_.d_min(d_min); }

			template <typename Vector>
			void set_d_max(Vector const& d_max) { qp_.d_max(d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { lowerBound_.x(x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { upperBound_.x(x_max); }

			template <typename Vector>
			void set_u_min(Vector const& u_min)
			{
				lowerBound_.u(u_min);
			}

			template <typename Vector>
			void set_u_min(size_t i, Vector const& u_min)
			{
				lowerBound_.u(i, u_min);
			}

			template <typename Vector>
			void set_u_max(Vector const& u_max)
			{
				upperBound_.u(u_max);
			}

			template <typename Vector>
			void set_u_max(size_t i, Vector const& u_max)
			{
				upperBound_.u(i, u_max);
			}

		private:
			std::size_t const i_;
			MpcRealtimeIteration& rti_;
			QpStageBase<typename QpWorkspace::Stage>& qp_;

			/// \brief The working point (linearization point).
			TrajectoryPoint& workingPoint_;

			/// \brief Lower bound of the working point.
			TrajectoryPoint& lowerBound_;

			/// \brief Upper bound of the working point.
			TrajectoryPoint& upperBound_;
		};


		/**
		 * \brief Set all data of the QP to NaN.
		 */
		void SetQpNaN()
		{
			for (auto& stage : qp_.problem())
			{
				stage.Q(nan());
				stage.R(nan());
				stage.S(nan());
				stage.q(nan());
				stage.r(nan());

				stage.A(nan());
				stage.B(nan());
				stage.b(nan());

				stage.lbx(nan());
				stage.ubx(nan());
				stage.lbu(nan());
				stage.ubu(nan());

				stage.C(nan());
				stage.D(nan());

				stage.lbd(nan());
				stage.ubd(nan());
			}
		}

		// Uses ocp_ to initialize _QP based on current working point workingPoint_.
		void InitQP()
		{
			SetQpNaN();

			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, *this);
				ocp_.InitStage(stage);
			}

			Stage terminal_stage(nT(), *this);
			ocp_.InitTerminalStage(terminal_stage);
		}

		// Uses ocp_ to update _QP based on current working point workingPoint_.
		void UpdateQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, *this);
				ocp_.UpdateStage(stage);
			}

			Stage terminal_stage(nT(), *this);
			ocp_.UpdateTerminalStage(terminal_stage);

#if (false && NDEBUG)
			// Error catching: check if we have any unexpected NaNs of infs in QP after update.

			std::vector<std::string> messages;
			qp::diagnose(qp_, std::back_inserter(messages));

			if (!messages.empty())
			{
				std::stringstream err_msg;
				err_msg << "MpcRealtimeIteration encountered an ill-formed QP. Messages follow:" << std::endl;
				std::copy(messages.begin(), messages.end(), std::ostream_iterator<std::string>(err_msg, "\n"));

				throw std::logic_error(err_msg.str());
			}
#endif
		}

		Real levenbergMarquardt_ = 0.;

		/// \brief The working point (linearization point).
		WorkingPoint workingPoint_;

		/// \brief Lower bound of the working point.
		WorkingPoint lowerBound_;

		/// \brief Upper bound of the working point.
		WorkingPoint upperBound_;

		/// \brief Multistage QP workspace.
		QpWorkspace qp_;

		// Private data members.
		OCP const& ocp_;	// <-- TODO: change to value semantics
		
		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool prepared_;

		static Real constexpr inf()
		{
			return std::numeric_limits<Real>::infinity();
		}

		static Real constexpr nan()
		{
			return std::numeric_limits<Real>::signaling_NaN();
		}
	};
}
