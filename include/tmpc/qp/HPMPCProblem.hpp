#pragma once

#include <tmpc/Matrix.hpp>
#include "QpSize.hpp"

#include <vector>
#include <limits>

namespace tmpc
{
	/* \brief Stores data for a multistage QP problem.
	 * The storage format is what is expected by the c_order_d_ip_ocp_hard_tv() function from HPMPC.
	 *
	 * TODO: convert the formulas to LaTeX/Doxygen format
	 *
	 *  The problem is stated as following:
	*
	*	min  sum_{ k = 0..nI } z_k'*H_k*z_k + g_k'*z_k
	*	s.t. x_{ k + 1 } = C_k * z_k + c_k				for k = 0..nI - 1
	*            dLow_k <= D_k * z_k <= dUpp_k			for k = 0..nI
	*            zMin_k <= z_k <= zMax_k                for k = 0..nI
	*
	*	where x_k is implicitly defined by z_k = [x_k  u_k] as the first nX variables of z_k
	*
	*	It holds
	*	z_k  \in R^nZ  for k = 0..nI - 1
	*   z_nI \in R*nX
	*
	*	nX < nZ
	*	nU = nZ - nX
	*
	*	TODO: parameterize HPMPCProblem by a kernel class.
	*/
	template <typename Scalar_>
	class HPMPCProblem
	{
	public:
		using Scalar = Scalar_;

		class Stage
		{
		public:
			static auto constexpr storageOrder = rowMajor;
			typedef DynamicMatrix<Scalar, storageOrder> Matrix;
			typedef DynamicVector<Scalar, columnVector> Vector;

			Stage(QpSize const& sz, size_t nx_next)
			:	size_(sz)
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

			const Matrix& get_A() const {
				return A_;
			}

			template <typename T>
			void set_A(const T& a) {
				full(A_) = a;
			}

			const Matrix& get_B() const {
				return B_;
			}

			template <typename T>
			void set_B(const T& b) {
				full(B_) = b;
			}

			Vector const& get_b() const {
				return b_;
			}

			template <typename T>
			void set_b(const T& b) {
				full(b_) = b;
			}

			const Matrix& get_C() const {
				return C_;
			}

			template <typename T>
			void set_C(const T& c) {
				full(C_) = c;
			}

			const Matrix& get_D() const {
				return D_;
			}

			template <typename T>
			void set_D(const T& d) {
				full(D_) = d;
			}

			const Vector& get_lbd() const {
				return lbd_;
			}

			template <typename T>
			void set_lbd(const T& lbd) {
				full(lbd_) = lbd;
			}

			Subvector<Vector const> get_lbu() const {
				return subvector(lb_, 0, size_.nu());
			}

			template <typename T>
			void set_lbu(const T& lbu) {
				subvector(lb_, 0, size_.nu()) = lbu;
			}

			Subvector<Vector const> get_lbx() const {
				return subvector(lb_, size_.nu(), size_.nx());
			}

			template <typename T>
			void set_lbx(const T& lbx) {
				subvector(lb_, size_.nu(), size_.nx()) = lbx;
			}

			const Matrix& get_Q() const {
				return Q_;
			}

			template <typename T>
			void set_Q(const T& q) {
				full(Q_) = q;
			}

			const Matrix& get_R() const {
				return R_;
			}

			template <typename T>
			void set_R(const T& r) {
				full(R_) = r;
			}

			const Matrix& get_S() const {
				return S_;
			}

			template <typename T>
			void set_S(const T& s) {
				full(S_) = s;
			}

			const Vector& get_q() const {
				return q_;
			}

			template <typename T>
			void set_q(const T& q) {
				full(q_) = q;
			}

			const Vector& get_r() const {
				return r_;
			}

			template <typename T>
			void set_r(const T& r) {
				full(r_) = r;
			}

			const Vector& get_ubd() const {
				return ubd_;
			}

			template <typename T>
			void set_ubd(const T& ubd) {
				full(ubd_) = ubd;
			}

			Subvector<Vector const> get_ubu() const {
				return subvector(ub_, 0, size_.nu());
			}

			template <typename T>
			void set_ubu(const T& ubu) {
				subvector(ub_, 0, size_.nu()) = ubu;
			}

			Subvector<Vector const> get_ubx() const {
				return subvector(ub_, size_.nu(), size_.nx());
			}

			template <typename T>
			void set_ubx(const T& ubx) {
				subvector(ub_, size_.nu(), size_.nx()) = ubx;
			}

			QpSize const& size() const
			{
				return size_;
			}

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
		};

		Stage& operator[](std::size_t i)
		{
			return stage_.at(i);
		}

		Stage const& operator[](std::size_t i) const
		{
			return stage_.at(i);
		}

		std::size_t size() const
		{
			return stage_.size();
		}

		typedef typename std::vector<Stage>::iterator iterator;
		typedef typename std::vector<Stage>::const_iterator const_iterator;
		typedef typename std::vector<Stage>::reference reference;
		typedef typename std::vector<Stage>::const_reference const_reference;

		iterator begin()
		{
			return stage_.begin();
		}

		iterator end()
		{
			return stage_.end();
		}

		const_iterator begin() const
		{
			return stage_.begin();
		}

		const_iterator end() const
		{
			return stage_.end();
		}

		reference front()
		{
			return stage_.front();
		}

		reference back()
		{
			return stage_.back();
		}

		const_reference front() const
		{
			return stage_.front();
		}

		const_reference back() const
		{
			return stage_.back();
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

		template <typename InputIterator>
		HPMPCProblem(InputIterator size_first, InputIterator size_last)
		{
			auto const nt = std::distance(size_first, size_last);
			stage_.reserve(nt);

			_nx.reserve(nt);
			_nu.reserve(nt);
			_nb.reserve(nt);
			_ng.reserve(nt);

			for (auto sz = size_first; sz != size_last; ++sz)
			{
				stage_.emplace_back(*sz, sz + 1 != size_last ? sz[1].nx() : 0);
				auto& st = stage_.back();

				_nx.push_back(sz->nx());
				_nu.push_back(sz->nu());
				_nb.push_back(sz->nx() + sz->nu());
				_ng.push_back(sz->nc());
			}

			_A .resize(nt);
			_B .resize(nt);
			_b .resize(nt);
			_Q .resize(nt);
			_S .resize(nt);
			_R .resize(nt);
			_q .resize(nt);
			_r .resize(nt);
			_lb.resize(nt);
			_ub.resize(nt);
			_C .resize(nt);
			_D .resize(nt);
			_lg.resize(nt);
			_ug.resize(nt);

			InitPointers();
		}

		HPMPCProblem(HPMPCProblem const& rhs)
		:	stage_(rhs.stage_)
		,	_A(rhs.size())
		,	_B(rhs.size())
		,	_b(rhs.size())
		,	_Q(rhs.size())
		,	_S(rhs.size())
		,	_R(rhs.size())
		,	_q(rhs.size())
		,	_r (rhs.size())
		,	_lb(rhs.size())
		,	_ub(rhs.size())
		,	_C (rhs.size())
		,	_D (rhs.size())
		,	_lg(rhs.size())
		,	_ug(rhs.size())
		,	_nx(rhs._nx)
		,	_nu(rhs._nu)
		,	_nb(rhs._nb)
		,	_ng(rhs._ng)
		{
			InitPointers();
		}

		HPMPCProblem(HPMPCProblem && rhs) = default;

	private:
		void InitPointers()
		{
			// Filling out pointer arrays for HPMPC
			for (size_t i = 0; i < size(); ++i)
			{
				_A[i] = stage_[i].A_data();
				_B[i] = stage_[i].B_data();
				_b[i] = stage_[i].b_data();

				_Q [i] = stage_[i].Q_data();
				_S [i] = stage_[i].S_data();
				_R [i] = stage_[i].R_data();
				_q [i] = stage_[i].q_data();
				_r [i] = stage_[i].r_data();
				_lb[i] = stage_[i].lb_data();
				_ub[i] = stage_[i].ub_data();

				_C [i] = stage_[i].C_data();
				_D [i] = stage_[i].D_data();
				_lg[i] = stage_[i].lg_data();
				_ug[i] = stage_[i].ug_data();
			}
		}

		// Private data members.
		//
		static Scalar constexpr sNaN()
		{
			return std::numeric_limits<Scalar>::signaling_NaN();
		}

		// Stores stage data
		std::vector<Stage> stage_;

		// "A" data array for HPMPC
		std::vector<double const *> _A;

		// "B" data array for HPMPC
		std::vector<double const *> _B;

		// "b" data array for HPMPC
		std::vector<double const *> _b;

		// "Q" data array for HPMPC
		std::vector<double const *> _Q;

		// "S" data array for HPMPC
		std::vector<double const *> _S;

		// "R" data array for HPMPC
		std::vector<double const *> _R;

		// "q" data array for HPMPC
		std::vector<double const *> _q;

		// "r" data array for HPMPC
		std::vector<double const *> _r;

		// "lb" data array for HPMPC
		std::vector<double const *> _lb;

		// "ub" data array for HPMPC
		std::vector<double const *> _ub;

		// "C" data array for HPMPC
		std::vector<double const *> _C;

		// "D" data array for HPMPC
		std::vector<double const *> _D;

		// "lg" data array for HPMPC
		std::vector<double const *> _lg;

		// "ug" data array for HPMPC
		std::vector<double const *> _ug;

		// Array of NX sizes
		std::vector<int> _nx;

		// Array of NU sizes
		std::vector<int> _nu;

		// Array of NB (bound constraints) sizes
		std::vector<int> _nb;

		// Array of NG (path constraints) sizes
		std::vector<int> _ng;
	};
}
