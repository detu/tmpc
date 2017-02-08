#pragma once

#include "../core/matrix.hpp"
#include "qp.hpp"
#include "QpSize.hpp"
#include "../core/RealtimeIteration.hpp"	// for RtiQpSize only

#include <Eigen/Dense>
#include <blaze/Math.h>

#include <vector>
#include <initializer_list>

namespace tmpc {

/** \brief Stores data for a multistage QP problem.
 *
 * Storage format is not explicitly defined and no access to raw data is provided..
 *
 *	TODO: Update problem statement.
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
*/
template <unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
class QuadraticProblemEigen
{
public:
	typedef unsigned int size_type;

	static unsigned const NX = NX_;
	static unsigned const NU = NU_;
	static unsigned const NC = NC_;
	static unsigned const NCT = NCT_;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;
	typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
	typedef Eigen::Matrix<double, NU, NU> InputInputMatrix;
	typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;
	typedef Eigen::Matrix<double, NC, NX> ConstraintStateMatrix;
	typedef Eigen::Matrix<double, NC, NU> ConstraintInputMatrix;
	typedef Eigen::Matrix<double, NC, 1> ConstraintVector;
	typedef Eigen::Matrix<double, NCT, NX> EndConstraintStateMatrix;
	typedef Eigen::Matrix<double, NCT, 1> EndConstraintVector;

	size_type nT() const { return _nt; }
	static constexpr size_type nX() { return NX; }
	static constexpr size_type nU() { return NU; }
	static constexpr size_type nD() { return NC; }
	static constexpr size_type nDT() { return NCT; }

	StateStateMatrix const& get_Q(size_type i) const { return stage(i, 1).Q; }
	template <typename Matrix> void set_Q(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).Q = val; }

	/**
	 * \brief Get R matrix of stage k
	 */
	InputInputMatrix const& get_R(size_type k) const
	{
		return stage(k).R;
	}

	/**
	 * \brief Set R matrix of a given stage
	 */
	template <typename Matrix>
	void set_R(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).R = val;
	}

	/**
	 * \brief Set block of R matrix of stage k with its top left corner at (i, j)
	 */
	template <typename Matrix>
	void set_R(size_type k, unsigned i, unsigned j, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).R.template block<Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime>(i, j) = val;
	}

	/**
	 * \brief Get S matrix of stage k
	 */
	StateInputMatrix const& get_S(size_type k) const
	{
		return stage(k).S;
	}

	/**
	 * \brief Set S matrix of stage k
	 */
	template <typename Matrix>
	void set_S(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).S = val;
	}

	/**
	 * \brief Set block of S matrix of stage k with its top left corner at (i, j)
	 */
	template <typename Matrix>
	void set_S(size_type k, unsigned i, unsigned j, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).S.template block<Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime>(i, j) = val;
	}

	/**
	*/
	StateVector const& get_q(size_type i) const { return stage(i, 1).q; }
	template <typename Matrix> void set_q(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).q = val; }

	/**
	 * \brief Get r vector of stage k
	 */
	InputVector const& get_r(size_type k) const
	{
		return stage(k).r;
	}

	/**
	 * \brief Set r vector of stage k
	 */
	template <typename Matrix>
	void set_r(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).r = val;
	}

	/**
	 * \brief Set a block of r vector of stage k starting from element i
	 */
	template <typename Matrix>
	void set_r(size_type k, unsigned i, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).r.template middleRows<Matrix::RowsAtCompileTime>(i) = val;
	}

	/**
	*/
	StateStateMatrix const& get_A(size_type i) const { return stage(i).A; }
	template <typename Matrix> void set_A(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).A = val; }

	/**
	 * \brief Get B matrix of stage k
	 */
	StateInputMatrix const& get_B(size_type k) const
	{
		return stage(k).B;
	}

	/**
	 * \brief Set B matrix of stage k
	 */
	template <typename Matrix>
	void set_B(size_type k, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).B = val;
	}

	/**
	 * \brief Set a block of B matrix of stage k with its top left corner at (i, j)
	 */
	template <typename Matrix>
	void set_B(size_type k, unsigned i, unsigned j, Eigen::MatrixBase<Matrix> const& val)
	{
		stage(k).B.template block<Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime>(i, j) = val;
	}

	/**
	*/
	StateVector const& get_b(size_type i) const { return stage(i).b; }
	template <typename Matrix> void set_b(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).b = val; }

	ConstraintStateMatrix const& get_C(size_type i) const { return stage(i).C; }
	template <typename Matrix> void set_C(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).C = val; }

	ConstraintInputMatrix const& get_D(size_type i) const { return stage(i).D; }
	template <typename Matrix> void set_D(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).D = val; }

	EndConstraintStateMatrix const& get_C_end() const { return _C_end; }
	template <typename Matrix> void set_C_end(Eigen::MatrixBase<Matrix> const& val) { _C_end = val; }

	ConstraintVector const& get_d_min(size_type i) const { return stage(i).d_min; }
	template <typename Matrix> void set_d_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).d_min = val; }

	EndConstraintVector const& get_d_end_min() const { return _d_end_min; }
	template <typename Matrix> void set_d_end_min(Eigen::MatrixBase<Matrix> const& val)	{ _d_end_min = val; }

	ConstraintVector const& get_d_max(size_type i) const { return stage(i).d_max; }
	template <typename Matrix> void set_d_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).d_max = val; }

	EndConstraintVector const& get_d_end_max() const { return _d_end_max; }
	template <typename Matrix> void set_d_end_max(Eigen::MatrixBase<Matrix> const& val)	{ _d_end_max = val; }

	StateVector const& get_x_min(size_type i) const { return stage(i, 1).x_min; }
	template <typename Matrix> void set_x_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).x_min = val; }

	StateVector const& get_x_max(size_type i) const { return stage(i, 1).x_max; }
	template <typename Matrix> void set_x_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i, 1).x_max = val; }

	InputVector const& get_u_min(size_type i) const { return stage(i).u_min; }
	template <typename Matrix> void set_u_min(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).u_min = val; }

	InputVector const& get_u_max(size_type i) const { return stage(i).u_max; }
	template <typename Matrix> void set_u_max(size_type i, Eigen::MatrixBase<Matrix> const& val) { stage(i).u_max = val; }

	std::vector<QpSize> const& size() const
	{
		return size_;
	}

	QuadraticProblemEigen(size_type nt)
	:	_nt(nt)
	,	size_(RtiQpSize(nt, NX, NU, NC, NCT))
	,	_stage(nt + 1)	// only some fields of _stage[nt] are used
	{}

	// Default copy constructor is ok.
	QuadraticProblemEigen(QuadraticProblemEigen const&) = default;

private:
	struct StageData
	{
		StateStateMatrix Q = signaling_nan<StateStateMatrix>();
		StateInputMatrix S = signaling_nan<StateInputMatrix>();
		InputInputMatrix R = signaling_nan<InputInputMatrix>();
		StateVector q = signaling_nan<StateVector>();
		InputVector r = signaling_nan<InputVector>();
		StateStateMatrix A = signaling_nan<StateStateMatrix>();
		StateInputMatrix B = signaling_nan<StateInputMatrix>();
		StateVector b = signaling_nan<StateVector>();
		StateVector x_min = signaling_nan<StateVector>();
		StateVector x_max = signaling_nan<StateVector>();
		InputVector u_min = signaling_nan<InputVector>();
		InputVector u_max = signaling_nan<InputVector>();
		ConstraintStateMatrix C = signaling_nan<ConstraintStateMatrix>();
		ConstraintInputMatrix D = signaling_nan<ConstraintInputMatrix>();
		ConstraintVector d_min = signaling_nan<ConstraintVector>();
		ConstraintVector d_max = signaling_nan<ConstraintVector>();
	};

	StageData& stage(size_type i, std::size_t delta = 0)
	{
		assert(i < nT() + delta && i < _stage.size());
		return _stage[i];
	}

	StageData const& stage(size_type i, std::size_t delta = 0) const
	{
		assert(i < nT() + delta && i < _stage.size());
		return _stage[i];
	}

	// Private data members.
	//

	// Number of control intervals.
	std::size_t _nt;

	std::vector<QpSize> size_;

	// Stores stage data
	std::vector<StageData> _stage;

	// 1 matrix of size _Nx x _Nx.
	StateStateMatrix _Q_end = signaling_nan<StateStateMatrix>();

	// 1 vector of size _Nx
	StateVector _q_end = signaling_nan<StateVector>();

	// 1 matrix of size NCT x NX.
	EndConstraintStateMatrix _C_end = signaling_nan<EndConstraintStateMatrix>();

	// 1 vector of size NCT
	EndConstraintVector _d_end_min = signaling_nan<EndConstraintVector>();

	// 1 vector of size NCT
	EndConstraintVector _d_end_max = signaling_nan<EndConstraintVector>();

	// 1 vector of size _Nx
	StateVector _x_end_min = signaling_nan<StateVector>();

	// 1 vector of size _Nx
	StateVector _x_end_max = signaling_nan<StateVector>();
};

/** \brief Stores data for a multistage QP problem.
 *
 * Storage format is not explicitly defined and no access to raw data is provided..
 *
 *	TODO: Update problem statement.
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
*/
template <typename Scalar_>
class QuadraticProblemBlaze
{
public:
	typedef std::size_t size_type;
	typedef Scalar_ Scalar;
	typedef blaze::DynamicMatrix<Scalar> Matrix;
	typedef blaze::DynamicVector<Scalar> Vector;

	class Stage
	{
	public:
		Stage(QpSize const& sz, size_type nx_next)
		:	size_(sz)
		,	Q_(sz.nx(), sz.nx())
		,	R_(sz.nu(), sz.nu())
		,	S_(sz.nx(), sz.nu())
		,	q_(sz.nx())
		,	r_(sz.nu())
		,	A_(nx_next, sz.nx())
		,	B_(nx_next, sz.nu())
		,	b_(nx_next)
		,	C_(sz.nc(), sz.nx())
		,	D_(sz.nc(), sz.nu())
		,	lbx_(sz.nx())
		,	ubx_(sz.nx())
		,	lbu_(sz.nu())
		,	ubu_(sz.nu())
		,	lbd_(sz.nc())
		,	ubd_(sz.nc())
		{
		}

		Stage(Stage const&) = delete;
		Stage(Stage &&) = default;

		template <typename Expr>
		Stage& operator=(Expr const& rhs)
		{
			assign(*this, rhs);
			return *this;
		}

		const Matrix& Q() const { return Q_; }
		Matrix& Q() { return Q_; }

		const Matrix& R() const { return R_; }
		Matrix& R() { return R_; }

		const Matrix& S() const { return S_; }
		Matrix& S() { return S_; }

		const Vector& q() const { return q_; }
		Vector& q() { return q_; }

		const Vector& r() const { return r_; }
		Vector& r() { return r_; }

		const Matrix& A() const { return A_; }
		Matrix& A() { return A_; }

		const Matrix& B() const { return B_; }
		Matrix& B() { return B_; }

		Vector const& b() const { return b_; }
		Vector& b() { return b_; }

		const Matrix& C() const { return C_; }
		Matrix& C() { return C_; }

		const Matrix& D() const { return D_; }
		Matrix& D() { return D_; }

		const Vector& lbx() const {	return lbx_; }
		Vector& lbx() {	return lbx_; }

		const Vector& ubx() const {	return ubx_; }
		Vector& ubx() {	return ubx_; }

		const Vector& lbu() const { return lbu_; }
		Vector& lbu() {	return lbu_; }

		const Vector& ubu() const {	return ubu_; }
		Vector& ubu() {	return ubu_; }

		const Vector& lbd() const {	return lbd_; }
		Vector& lbd() {	return lbd_; }

		const Vector& ubd() const {	return ubd_; }
		Vector& ubd() {	return ubd_; }

		QpSize const& size() const { return size_; }

	private:
		QpSize size_;
		Matrix Q_;
		Matrix R_;
		Matrix S_;
		Vector q_;
		Vector r_;
		Matrix A_;
		Matrix B_;
		Vector b_;
		Matrix C_;
		Matrix D_;
		Vector lbx_;
		Vector ubx_;
		Vector lbu_;
		Vector ubu_;
		Vector lbd_;
		Vector ubd_;
	};

	Stage& operator[](std::size_t i)
	{
		return stage_.at(i);
	}

	Stage const& operator[](std::size_t i) const
	{
		return stage_.at(i);
	}

	Stage& stage(std::size_t i)
	{
		return stage_.at(i);
	}

	Stage const& stage(std::size_t i) const
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

	QuadraticProblemBlaze(std::initializer_list<QpSize> sz)
	:	QuadraticProblemBlaze(sz.begin(), sz.end())
	{
	}

	template <typename InIter>
	QuadraticProblemBlaze(InIter sz_begin, InIter sz_end)
	{
		stage_.reserve(std::distance(sz_begin, sz_end));

		for (auto sz = sz_begin; sz != sz_end; ++sz)
		{
			auto const nx_next = sz + 1 != sz_end ? (sz + 1)->nx() : 0;
			stage_.emplace_back(*sz, nx_next);
		}
	}

	// Default copy constructor is ok.
	QuadraticProblemBlaze(QuadraticProblemBlaze const&) = default;

private:
	// Private data members.
	//

	// Stores stage data
	std::vector<Stage> stage_;
};

}	// namespace tmpc
