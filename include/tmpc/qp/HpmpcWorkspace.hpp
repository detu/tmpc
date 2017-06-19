#pragma once

#include "QpSolverException.hpp"
#include "QpSize.hpp"

#include <tmpc/Matrix.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>

namespace tmpc
{
	namespace hpmpc_wrapper
	{
		int c_order_d_ip_ocp_hard_tv(	int *kk, int k_max, double mu0, double mu_tol,
										int N, int const *nx, int const *nu, int const *nb, int const *ng,
										int warm_start,
										double const * const *A, double const * const *B, double const * const *b,
										double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
										double const * const *lb, double const * const *ub,
										double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
										double * const *x, double * const *u, double * const *pi, double * const *lam, double * const *t,
										double *inf_norm_res,
										void *work0,
										double *stat );

		int d_ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const *ng);
	}

	class HpmpcException : public QpSolverException
	{
	public:
		HpmpcException(int code)
		:	QpSolverException("HPMPC"),
			_code(code),
			msg_(std::string(QpSolverException::what()) + "\nReturn code " + std::to_string(code))
		{
		}

		int getCode() const	{ return _code;	}
		char const * what() const noexcept override { return msg_.c_str(); }

	private:
		int const _code;
		std::string const msg_;
	};

	/**
	 * \brief Multistage QP solver using qpOASES
	 *
	 * \tparam <Scalar_> Scalar type
	 */
	template <typename Scalar_>
	class HpmpcWorkspace
	{
	public:
		typedef Scalar_ Scalar;

		class Stage
		{
		public:
			static auto constexpr storageOrder = rowMajor;
			typedef DynamicMatrix<Scalar, storageOrder> Matrix;
			typedef DynamicVector<Scalar, columnVector> Vector;

			Stage(QpSize const& sz, size_t nx_next)
			:	size_(sz)
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			,	Q_(sz.nx(), sz.nx(), sNaN())
			,	R_(sz.nu(), sz.nu(), sNaN())
			,	S_(sz.nx(), sz.nu(), sNaN())
			,	q_(sz.nx(), sNaN())
			,	r_(sz.nu(), sNaN())
			,	A_(nx_next, sz.nx(), sNaN())
			,	B_(nx_next, sz.nu(), sNaN())
			,	b_(nx_next, sNaN())
			,	C_(sz.nc(), sz.nx(), sNaN())
			,	D_(sz.nc(), sz.nu(), sNaN())
			,	lb_(sz.nu() + sz.nx(), sNaN())
			,	ub_(sz.nu() + sz.nx(), sNaN())
			,	lbd_(sz.nc(), sNaN())
			,	ubd_(sz.nc(), sNaN())
			,	_x(sz.nx(), sNaN())
			,	_u(sz.nu(), sNaN())
			,	_pi(nx_next, sNaN())
			,	_lam(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN())
			,	_t(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN())
			{
			}

			Stage(Stage const&) = default; //delete;
			Stage(Stage &&) = default;

			template <typename Expr>
			Stage& operator=(Expr const& rhs)
			{
				assign(*this, rhs);
				return *this;
			}

			const Matrix& A() const { return A_; }
			template <typename T> void A(const T& a) { full(A_) = a; }

			const Matrix& B() const { return B_; }
			template <typename T> void B(const T& b) { full(B_) = b; }

			Vector const& b() const { return b_; }
			template <typename T> void b(const T& b) { full(b_) = b; }

			const Matrix& C() const { return C_; }
			template <typename T> void C(const T& c) { full(C_) = c; }

			const Matrix& D() const { return D_; }
			template <typename T> void D(const T& d) { full(D_) = d; }

			const Vector& lbd() const {	return lbd_; }
			template <typename T> void lbd(const T& lbd) { full(lbd_) = lbd; }

			Subvector<Vector const> lbu() const { return subvector(lb_, 0, size_.nu());	}			
			template <typename T> void lbu(const T& lbu) { subvector(lb_, 0, size_.nu()) = lbu; }

			Subvector<Vector const> lbx() const { return subvector(lb_, size_.nu(), size_.nx()); }
			template <typename T> void lbx(const T& lbx) { subvector(lb_, size_.nu(), size_.nx()) = lbx; }

			const Matrix& Q() const { return Q_; }
			template <typename T> void Q(const T& q) { full(Q_) = q; }

			const Matrix& R() const { return R_; }
			template <typename T> void R(const T& r) { full(R_) = r; }

			const Matrix& S() const { return S_; }
			template <typename T> void S(const T& s) { full(S_) = s; }

			const Vector& q() const { return q_; }
			template <typename T> void q(const T& q) { full(q_) = q; }

			const Vector& r() const { return r_; }
			template <typename T> void r(const T& r) { full(r_) = r; }

			const Vector& ubd() const { return ubd_; }
			template <typename T> void ubd(const T& ubd) { full(ubd_) = ubd; }

			Subvector<Vector const> ubu() const { return subvector(ub_, 0, size_.nu()); }
			template <typename T> void ubu(const T& ubu) { subvector(ub_, 0, size_.nu()) = ubu; }

			Subvector<Vector const> ubx() const { return subvector(ub_, size_.nu(), size_.nx()); }
			template <typename T> void ubx(const T& ubx) { subvector(ub_, size_.nu(), size_.nx()) = ubx; }

			DynamicVector<Scalar> const& x() const { return _x; }
			DynamicVector<Scalar> const& u() const { return _u;	}
			DynamicVector<Scalar> const& pi() const	{ return _pi; }
			DynamicVector<Scalar> const& lam() const { return _lam; }
			DynamicVector<Scalar> const& t() const { return _t; }

			QpSize const& size() const { return size_; }

			// ******************************************************
			//                HPMPC raw data interface.
			//
			// The prefixes before _data() correspond to the names of
			// the argument to c_order_d_ip_ocp_hard_tv().
			// ******************************************************
			double const * A_data () const { return A_.data(); }
			double const * B_data () const { return B_.data(); }
			double const * b_data () const { return b_.data();	}
			double const * Q_data () const { return Q_.data(); }
			double const * S_data () const { return S_.data(); }
			double const * R_data () const { return R_.data(); }
			double const * q_data () const { return q_.data(); }
			double const * r_data () const { return r_.data();	}
			double const * lb_data() const { return lb_.data(); }
			double const * ub_data() const { return ub_.data(); }
			double const * C_data () const { return C_.data(); }
			double const * D_data () const { return D_.data(); }
			double const * lg_data() const { return lbd_.data(); }
			double const * ug_data() const { return ubd_.data(); }

			double * x_data() { return _x.data(); }
			double * u_data() { return _u.data(); }
			double * pi_data() { return _pi.data(); }
			double * lam_data() { return _lam.data(); }
			double * t_data() { return _t.data(); }

		private:
			QpSize size_;

			// Hessian = [R, S; S', Q]
			Matrix R_;
			Matrix S_;
			Matrix Q_;

			// Gradient = [r; q]
			Vector r_;
			Vector q_;

			// Inter-stage equalities x_{k+1} = A x_k + B u_k + c_k
			Matrix A_;
			Matrix B_;
			Vector b_;

			// Inequality constraints d_{min} <= C x_k + D u_k <= d_{max}
			Matrix C_;
			Matrix D_;
			Vector lbd_;
			Vector ubd_;

			// Bound constraints:
			// lb <= [u; x] <= ub
			Vector lb_;
			Vector ub_;

			DynamicVector<Scalar> _x;
			DynamicVector<Scalar> _u;
			DynamicVector<Scalar> _pi;
			DynamicVector<Scalar> _lam;
			DynamicVector<Scalar> _t;
		};

		Stage& operator[](std::size_t i) { return stage_.at(i); }
		Stage const& operator[](std::size_t i) const { return stage_.at(i);	}

		std::size_t size() const { return stage_.size(); }

		typedef typename std::vector<Stage>::iterator iterator;
		typedef typename std::vector<Stage>::const_iterator const_iterator;
		typedef typename std::vector<Stage>::reference reference;
		typedef typename std::vector<Stage>::const_reference const_reference;

		iterator begin() { return stage_.begin(); }
		iterator end() { return stage_.end(); }

		const_iterator begin() const { return stage_.begin(); }
		const_iterator end() const { return stage_.end(); }

		reference front() { return stage_.front(); }
		reference back() { return stage_.back(); }

		const_reference front() const { return stage_.front(); }
		const_reference back() const { return stage_.back(); }

		HpmpcWorkspace(int max_iter = 100)
		:	_stat(max_iter)
		{
		}

		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename InputIterator>
		HpmpcWorkspace(InputIterator size_first, InputIterator size_last, int max_iter = 100)
		:	_stat(max_iter)
		,	_inf_norm_res(sNaN())
		{
			auto const nt = std::distance(size_first, size_last);
			stage_.reserve(nt);

			_nx.reserve(nt);
			_nu.reserve(nt);
			_nb.reserve(nt);
			_ng.reserve(nt);

			_A .reserve(nt);
			_B .reserve(nt);
			_b .reserve(nt);
			_Q .reserve(nt);
			_S .reserve(nt);
			_R .reserve(nt);
			_q .reserve(nt);
			_r .reserve(nt);
			_lb.reserve(nt);
			_ub.reserve(nt);
			_C .reserve(nt);
			_D .reserve(nt);
			_lg.reserve(nt);
			_ug.reserve(nt);
			
			_x.reserve(nt);
			_u.reserve(nt);
			_pi.reserve(nt);
			_lam.reserve(nt);
			_t.reserve(nt);

			for (auto sz = size_first; sz != size_last; ++sz)
			{
				stage_.emplace_back(*sz, sz + 1 != size_last ? sz[1].nx() : 0);
				auto& st = stage_.back();

				_nx.push_back(sz->nx());
				_nu.push_back(sz->nu());
				_nb.push_back(sz->nx() + sz->nu());
				_ng.push_back(sz->nc());

				_A.push_back(st.A_data());
				_B.push_back(st.B_data());
				_b.push_back(st.b_data());

				_Q.push_back(st.Q_data());
				_S.push_back(st.S_data());
				_R.push_back(st.R_data());
				_q.push_back(st.q_data());
				_r.push_back(st.r_data());
				_lb.push_back(st.lb_data());
				_ub.push_back(st.ub_data());

				_C .push_back(st.C_data());
				_D .push_back(st.D_data());
				_lg.push_back(st.lg_data());
				_ug.push_back(st.ug_data());

				_x.push_back(st.x_data());
				_u.push_back(st.u_data());
				_pi.push_back(st.pi_data());
				_lam.push_back(st.lam_data());
				_t.push_back(st.t_data());
			}

			// TODO: calculate workspace size and call
			// workspace_.reserve();
		}

		/**
		 * \brief Copy constructor
		 *
		 * Copying is not allowed.
		 */
		HpmpcWorkspace(HpmpcWorkspace const&) = delete;

		/**
		 * \brief Move constructor
		 *
		 * Move-construction is ok.
		 */
		HpmpcWorkspace(HpmpcWorkspace&& rhs) = default;

		HpmpcWorkspace& operator= (HpmpcWorkspace const&) = delete;
		HpmpcWorkspace& operator= (HpmpcWorkspace &&) = delete;

		void solve()
		{
			if (size() > 0)
			{
				// Number of QP steps for HPMPC
				auto const N = size() - 1;

				// Make sure we have enough workspace.
				_workspace.resize(hpmpc_wrapper::d_ip_ocp_hard_tv_work_space_size_bytes(
						static_cast<int>(N), nx_data(), nu_data(), nb_data(), ng_data()));

				// Call HPMPC
				auto const ret = hpmpc_wrapper::c_order_d_ip_ocp_hard_tv(&numIter_, getMaxIter(), _mu0, _muTol, N,
						nx_data(), nu_data(), nb_data(), ng_data(), _warmStart ? 1 : 0, A_data(), B_data(), b_data(),
						Q_data(), S_data(), R_data(), q_data(), r_data(), lb_data(), ub_data(), C_data(), D_data(),
						lg_data(), ug_data(), x_data(), u_data(), pi_data(), lam_data(), t_data(), inf_norm_res_data(),
						_workspace.data(), _stat[0].data());

				if (ret != 0)
				{
					throw HpmpcException(ret);
				}
			}
		}

		std::size_t getMaxIter() const noexcept { return _stat.size(); }

		double getMuTol() const noexcept { return _muTol; }

		void setMuTol(double val)
		{
			if (val <= 0.)
				throw std::invalid_argument("mu tolerance for hpmpc must be positive");

			_muTol = val;
		}

		/// \brief Get number of iterations performed by the QP solver.
		unsigned getNumIter() const { return numIter_; }
		
	private:
		// --------------------------------
		//
		// HPMPC QP problem data
		//
		// --------------------------------

		// Stores stage data
		std::vector<Stage> stage_;

		// "A" data array for HPMPC
		std::vector<Scalar const *> _A;

		// "B" data array for HPMPC
		std::vector<Scalar const *> _B;

		// "b" data array for HPMPC
		std::vector<Scalar const *> _b;

		// "Q" data array for HPMPC
		std::vector<Scalar const *> _Q;

		// "S" data array for HPMPC
		std::vector<Scalar const *> _S;

		// "R" data array for HPMPC
		std::vector<Scalar const *> _R;

		// "q" data array for HPMPC
		std::vector<Scalar const *> _q;

		// "r" data array for HPMPC
		std::vector<Scalar const *> _r;

		// "lb" data array for HPMPC
		std::vector<Scalar const *> _lb;

		// "ub" data array for HPMPC
		std::vector<Scalar const *> _ub;

		// "C" data array for HPMPC
		std::vector<Scalar const *> _C;

		// "D" data array for HPMPC
		std::vector<Scalar const *> _D;

		// "lg" data array for HPMPC
		std::vector<Scalar const *> _lg;

		// "ug" data array for HPMPC
		std::vector<Scalar const *> _ug;

		// Array of NX sizes
		std::vector<int> _nx;

		// Array of NU sizes
		std::vector<int> _nu;

		// Array of NB (bound constraints) sizes
		std::vector<int> _nb;

		// Array of NG (path constraints) sizes
		std::vector<int> _ng;

		// --------------------------------
		//
		// HPMPC QP solution data
		//
		// --------------------------------

		std::vector<double *> _x;
		std::vector<double *> _u;
		std::vector<double *> _pi;
		std::vector<double *> _lam;
		std::vector<double *> _t;
		StaticVector<double, 4> _inf_norm_res;

		/// \brief Number of iterations performed by the QP solver.
		int numIter_ = 0;

		// --------------------------------
		//
		// HPMPC solver data
		//
		// --------------------------------

		// Workspace for HPMPC functions
		std::vector<char> _workspace;

		// Iteration statistics. HPMPC returns 5 double numbers per iteration.
		typedef std::array<double, 5> IterStat;
		std::vector<IterStat> _stat;

		double _mu0 = 0.;
		double _muTol = 1e-10;

		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;

		static Scalar constexpr sNaN()
		{
			return std::numeric_limits<Scalar>::signaling_NaN();
		}

		// ******************************************************
		//                HPMPC raw data interface.
		//
		// The prefixes before _data() correspond to the names of
		// the argument to c_order_d_ip_ocp_hard_tv().
		// ******************************************************
		double const * const * A_data () const { return _A .data(); }
		double const * const * B_data () const { return _B .data(); }
		double const * const * b_data () const { return _b.data();	}
		double const * const * Q_data () const { return _Q .data(); }
		double const * const * S_data () const { return _S .data(); }
		double const * const * R_data () const { return _R .data(); }
		double const * const * q_data () const { return _q .data(); }
		double const * const * r_data () const { return _r .data();	}
		double const * const * lb_data() const { return _lb.data(); }
		double const * const * ub_data() const { return _ub.data(); }
		double const * const * C_data () const { return _C .data(); }
		double const * const * D_data () const { return _D .data(); }
		double const * const * lg_data() const { return _lg.data(); }	// TODO: beware of x[0] being equality constrained!
		double const * const * ug_data() const { return _ug.data(); }	// TODO: beware of x[0] being equality constrained!
		int const * nx_data() const { return _nx.data(); }
		int const * nu_data() const { return _nu.data(); }
		int const * nb_data() const { return _nb.data(); }
		int const * ng_data() const { return _ng.data(); }

		double * const * x_data() { return _x.data(); }
		double * const * u_data() { return _u.data(); }
		double * const * pi_data() { return _pi.data(); }
		double * const * lam_data() { return _lam.data(); }
		double * const * t_data() { return _t.data(); }
		double * inf_norm_res_data() { return _inf_norm_res.data(); }
	};
}
