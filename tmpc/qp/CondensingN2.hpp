#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/Matrix.hpp>

#include <vector>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <array>

namespace tmpc
{
	/**
	 * \brief CondensingN2 algorithm with O(N^2) runtime in horizon length.
	 * 
	 * Implements the condensing algorithm with quadratic runtime as described in Section 5.4 of Joel's thesis: 
	 * ftp://ftp.esat.kuleuven.be/pub/stadius/ida/reports/13-217.pdf.
	 * 
	 * Manages resources needed for the condensing algorithm of a given size.
	 *
	 * \tparam <Kernel> Kernel for matrix arithmetic.
	 * 
	 * TODO: account for soft constraints.
	 */
	template <typename Kernel_>
	class CondensingN2
	{
	public:
		using Kernel = Kernel_;

		template <typename InIter>
		class Expression
		:	public OcpQpExpressionBase<Expression<InIter>>
		{
			using size_type = typename Kernel::size_t;
			using Matrix = DynamicMatrix<Kernel>;
			using Vector = DynamicVector<Kernel>;
			using Real = typename Kernel::Real;

		public:
			Expression(Expression const&) = delete;
			Expression(Expression &&) = default;

			Expression(CondensingN2& c, InIter first, InIter last)
			:	c_(c)
			,	first_ {first}
			,	last_ {last}
			{
				assert(first_ < last_);
			}


			template <typename QP>
			void implEvalTo(OcpQpBase<QP>& qp) const
			{
				auto const stages = boost::make_iterator_range(first_, last_);
				auto const N = stages.size();
				
				// Initialization of QP variables
				auto nu = first_->size().nu();
				auto nc = first_->size().nc();

				//K::left_cols(Sc_, nu) = first_->S();
				submatrix(c_.Sc_, 0, 0, c_.Sc_.rows(), nu) = first_->S();

				//K::top_rows(Cc_, nc) = first_->C();
				submatrix(c_.Cc_, 0, 0, nc, columns(c_.Cc_)) = first_->C();
				//K::top_left_corner(Dc_, nc, nu) = first_->D();
				submatrix(c_.Dc_, 0, 0, nc, nu) = first_->D();
				//K::head(lbd, nc) = first_->lbd();
				subvector(c_.lbd_, 0, nc) = first_->lbd();
				//K::head(ubd, nc) = first_->ubd();
				subvector(c_.ubd_, 0, nc) = first_->ubd();

				//K::head(lbu, nu) = first_->lbu();
				subvector(c_.lbu_, 0, nu) = first_->lbu();
				//K::head(ubu, nu) = first_->ubu();
				subvector(c_.ubu_, 0, nu) = first_->ubu();

				c_.A_ = first_->A();
				c_.B_ = first_->B();

				std::vector<DynamicVector<Kernel>> g(N + 1);				
				g[0] = DynamicVector<Kernel>(stages[0].size().nx(), 0);

				std::vector<DynamicMatrix<Kernel>> F(N + 1);

				// Precalculate F. Can it be optimized/eliminated?
				F[0] = IdentityMatrix<Kernel>(stages[0].size().nx());
				for (size_t i = 0; i < N; ++i)
					F[i + 1] = stages[i].A() * F[i];
	
				c_.Qc_ = 0;
				c_.qc_ = 0;
				//c_.Sc_ = 0;

				for (size_t i = 0; i < N; ++i)
				{
					// Update g
					g[i + 1] = stages[i].A() * g[i] + stages[i].b();
					
					// Update G
					std::vector<DynamicMatrix<Kernel>> G(N + 1);
					G[i + 1] = stages[i].B();
					
					for (size_t k = i + 1; k < N; ++k)
						G[k + 1] = stages[k].A() * G[k];	// <-- this can be sped-up by precomputing partial products of A_k

					// Update R
					DynamicMatrix<Kernel> W(rows(stages[N - 1].B()), columns(G[N]), 0);

					for (size_t k = N; k-- > i + 1; )
					{
						c_.R(k, i) = trans(stages[k].S()) * G[k] + trans(stages[k].B()) * W;
						c_.R(i, k) = trans(c_.R(k, i));
						W = stages[k].Q() * G[k] + trans(stages[k].A()) * W;
					}

					c_.R(i, i) = stages[i].R() + trans(stages[i].B()) * W;

					// Update Q
					c_.Qc_ += trans(F[i]) * stages[i].Q() * F[i];

					// Update q
					c_.qc_ += trans(F[i]) * (stages[i].Q() * g[i] + stages[i].q());

					// Update S
					c_.S(i) = trans(F[i]) * stages[i].S(); 
					for (size_t k = i + 1; k < N; ++k)
						c_.S(i) += trans(F[k]) * stages[k].Q() * G[k];
				}

				auto stage = stages.begin() + 1;
				for (size_t k = 1; k < N; ++k, ++stage)
				{
					auto const& sz = stage->size();
					// TODO: add nx1() to OcpSize
					auto const nx_next = stage->B().rows();

					// Update C
					//K::middle_rows(c_.Cc_, nc, sz.nx() + sz.nc()) << A, stage->C() * A;
					submatrix(c_.Cc_, nc          , 0, sz.nx(), columns(c_.Cc_)) = c_.A_;
					submatrix(c_.Cc_, nc + sz.nx(), 0, sz.nc(), columns(c_.Cc_)) = stage->C() * c_.A_;

					// Update D (which is initialized to 0)
					//K::middle_rows(c_.Dc_, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs_.nu() - nu),
					//									  stage->C() * B, stage->D(), K::zero(sz.nc(), cs_.nu() - nu - sz.nu());
					submatrix(c_.Dc_, nc          ,  0, sz.nx(),      nu) = c_.B_;
					submatrix(c_.Dc_, nc + sz.nx(),  0, sz.nc(),      nu) = stage->C() * c_.B_;
					submatrix(c_.Dc_, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->D();

					// Update lbd
					subvector(c_.lbd_, nc          , sz.nx()) = stage->lbx() - g[k];
					subvector(c_.lbd_, nc + sz.nx(), sz.nc()) = stage->lbd() - stage->C() * g[k];

					// Update ubd
					subvector(c_.ubd_, nc          , sz.nx()) = stage->ubx() - g[k];
					subvector(c_.ubd_, nc + sz.nx(), sz.nc()) = stage->ubd() - stage->C() * g[k];

					// Uplate lbu, ubu
					subvector(c_.lbu_, nu, sz.nu()) = stage->lbu();
					subvector(c_.ubu_, nu, sz.nu()) = stage->ubu();

					// Update A
					c_.A_ = stage->A() * c_.A_;

					// Update B
					{
						DynamicMatrix<Kernel> B_next(nx_next, nu + sz.nu());
						submatrix(B_next, 0,  0, nx_next,      nu) = stage->A() * c_.B_;
						submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->B();
						c_.B_ = std::move(B_next);
					}

					// Update indices
					nu += sz.nu();
					nc += sz.nx() + sz.nc();
				}

				// Calculate r by backward substitution
				{
					c_.r(N - 1) = stages[N - 1].r() + trans(stages[N - 1].S()) * g[N - 1];
					DynamicVector<Kernel> w = stages[N - 1].q() + stages[N - 1].Q() * g[N - 1];

					for (size_t k = N - 1; k --> 0 ;)
					{
						// Update r
						c_.r(k) = stages[k].r() + trans(stages[k].S()) * g[k] + trans(stages[k].B()) * w;
						w = stages[k].q() + stages[k].Q() * g[k] + trans(stages[k].A()) * w;
					}
				}

				// Set the state bounds of the condensed stage<tmpc/ocp/OcpSize.hpp>
				qp.lbx(first_->lbx());
				qp.ubx(first_->ubx());

				// Assign the output values.				
				qp.A(c_.A_);
				qp.B(c_.B_);
				qp.b(g[N]);				
				qp.C(c_.Cc_);
				qp.D(c_.Dc_);
				qp.lbd(c_.lbd_);
				qp.ubd(c_.ubd_);
				qp.lbu(c_.lbu_);
				qp.ubu(c_.ubu_);
				qp.Q(c_.Qc_);
				qp.R(c_.Rc_);
				qp.S(c_.Sc_);
				qp.q(c_.qc_);
				qp.r(c_.rc_);

				// TODO: implement soft constraints matrices recalculation.
				qp.Zl(sNaN<Real>());
				qp.Zu(sNaN<Real>());
				qp.zl(sNaN<Real>());
				qp.zu(sNaN<Real>());

				// TODO: what happens to the soft constraints index?
			}


			auto const& implSize() const
			{
				return c_.cs_;
			}


		private:
			CondensingN2& c_;
			InIter first_;
			InIter last_;
		};

		template <typename InIter>
		CondensingN2(InIter const& sz_first, InIter const& sz_last)
		:	cs_(condensedOcpSize(sz_first, sz_last))
		,	Qc_(cs_.nx(), cs_.nx())
		,   Rc_(cs_.nu(), cs_.nu())
		,   Sc_(cs_.nx(), cs_.nu())
		,   qc_(cs_.nx())
		,   rc_(cs_.nu())
		,   Cc_(cs_.nc(), cs_.nx())
		,   Dc_(cs_.nc(), cs_.nu(), 0)
		,   lbd_(cs_.nc())
		,   ubd_(cs_.nc())
		,   lbu_(cs_.nu())
		,   ubu_(cs_.nu())
		,	size_(sz_first, sz_last)
		{
			size_t const N = std::distance(sz_first, sz_last);
			//g_.resize(N);
			
			// Init cumulative sizes.
			cumNu_.reserve(N);
			//cumNc_.reserve(N);

			auto nu = 0;
			//auto nc = size_.front().nc();

			for (auto sz = size_.begin(); sz != size_.end(); ++sz)
			{
				cumNu_.push_back(nu);
				//cumNc_.push_back(nc);
				nu += sz->nu();
				//nc += sz->nx() + sz->nc();
			}
		}

		/**
		* \brief Perform the condensing.
		*/
		template <typename InIter>
		Expression<InIter> operator()(InIter first, InIter last)
		{
			if (condensedOcpSize(ocpSizeIterator(first), ocpSizeIterator(last)) != cs_)
				throw std::invalid_argument("QP size does not match the condensed size");

			return Expression<InIter> {*this, first, last};
		}


		auto const& condensedSize() const
		{
			return cs_;
		}
		

	private:
		OcpSize const cs_;

		DynamicMatrix<Kernel> Qc_;
		DynamicMatrix<Kernel> Rc_;
		DynamicMatrix<Kernel> Sc_;
		DynamicVector<Kernel> qc_;
		DynamicVector<Kernel> rc_;

		DynamicMatrix<Kernel> Cc_;
		DynamicMatrix<Kernel> Dc_;
		DynamicVector<Kernel> lbd_;
		DynamicVector<Kernel> ubd_;

		DynamicVector<Kernel> lbu_;
		DynamicVector<Kernel> ubu_;

		DynamicMatrix<Kernel> A_;
		DynamicMatrix<Kernel> B_;

		/*
		class PerStageData
		{
		public:
			PerStageData(CondensingN2& c, size_t nx, size_t nu, size_t nc, size_t ns)
			:	

		private:
			CondensingN2& 

			OcpSize stageSize_;
			OcpSize cumSize_;

			DynamicVector<Kernel> g_;
			Subvector<Kernel, DynamicVector<Kernel>> r_;
		};
		*/

		//std::vector<DynamicVector<Kernel>> g_;
		std::vector<OcpSize> size_;
		std::vector<size_t> cumNu_;
		//std::vector<size_t> cumNc_;


		decltype(auto) r(size_t k)
		{
			return subvector(rc_, cumNu_[k], size_[k].nu());
		}


		decltype(auto) R(size_t i, size_t j)
		{
			return submatrix(Rc_, cumNu_[i], cumNu_[j], size_[i].nu(), size_[j].nu());
		}


		decltype(auto) S(size_t j)
		{
			return submatrix(Sc_, 0, cumNu_[j], rows(Sc_), size_[j].nu());
		}
	};
}
