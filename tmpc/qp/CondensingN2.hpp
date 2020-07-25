#pragma once

#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSolution.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Exception.hpp>

#include <vector>
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
	 * \tparam <Real> Scalar floating point type.
	 * 
	 * TODO: account for soft constraints.
	 */
	template <typename Real_>
	class CondensingN2
	{
	public:
		using Real = Real_;

		template <OcpQp Qp, OcpSolution Solution>
		void operator()(Qp const& qp, Solution& sol) const
		{
			auto const stages = make_iterator_range(first_, last_);
			auto const N = stages.size();
			
			// Initialization of QP variables
			auto nu = first_->size().nu();
			auto nc = first_->size().nc();

			//K::left_cols(Sc_, nu) = first_->S();
			submatrix(Sc_, 0, 0, Sc_.rows(), nu) = first_->S();

			//K::top_rows(Cc_, nc) = first_->C();
			submatrix(Cc_, 0, 0, nc, columns(Cc_)) = first_->C();
			//K::top_left_corner(Dc_, nc, nu) = first_->D();
			submatrix(Dc_, 0, 0, nc, nu) = first_->D();
			//K::head(lbd, nc) = first_->lbd();
			subvector(lbd_, 0, nc) = first_->lbd();
			//K::head(ubd, nc) = first_->ubd();
			subvector(ubd_, 0, nc) = first_->ubd();

			//K::head(lbu, nu) = first_->lbu();
			subvector(lbu_, 0, nu) = first_->lbu();
			//K::head(ubu, nu) = first_->ubu();
			subvector(ubu_, 0, nu) = first_->ubu();

			B_ = first_->B();

			std::vector<DynamicVector<Real>> g(N + 1);				
			g[0] = DynamicVector<Real>(stages[0].size().nx(), 0);

			std::vector<DynamicMatrix<Real>> F(N + 1);

			// Precalculate F. Can it be optimized/eliminated?
			F[0] = IdentityMatrix<Kernel>(stages[0].size().nx());
			for (size_t i = 0; i < N; ++i)
				F[i + 1] = stages[i].A() * F[i];

			Qc_ = 0;
			qc_ = 0;
			//Sc_ = 0;

			for (size_t i = 0; i < N; ++i)
			{
				// Update g
				g[i + 1] = stages[i].A() * g[i] + stages[i].b();
				
				// Update G
				blockG_(i + 1, i) = stages[i].B();
				
				for (size_t k = i + 1; k < N; ++k)
					blockG_(k + 1, i) = stages[k].A() * blockG_(k, i);	// <-- this can be sped-up by precomputing partial products of A_k

				// Update R
				DynamicMatrix<Real> W(rows(stages[N - 1].B()), columns(blockG_(N, i)), 0); // <-- blockG_(N, ...) is undefined!

				for (size_t k = N; k-- > i + 1; )
				{
					R(k, i) = trans(stages[k].S()) * blockG_(k, i) + trans(stages[k].B()) * W;
					R(i, k) = trans(R(k, i));
					W = stages[k].Q() * blockG_(k, i) + trans(stages[k].A()) * W;
				}

				R(i, i) = stages[i].R() + trans(stages[i].B()) * W;

				// Update Q
				Qc_ += trans(F[i]) * stages[i].Q() * F[i];

				// Update q
				qc_ += trans(F[i]) * (stages[i].Q() * g[i] + stages[i].q());

				// Update S
				S(i) = trans(F[i]) * stages[i].S(); 
				for (size_t k = i + 1; k < N; ++k)
					S(i) += trans(F[k]) * stages[k].Q() * blockG_(k, i);

				// Update C
				if (i > 0)
					Cx(i) = F[i];
				Cc(i) = stages[i].C() * F[i];
			}

			auto stage = stages.begin() + 1;
			for (size_t k = 1; k < N; ++k, ++stage)
			{
				auto const& sz = stage->size();
				// TODO: add nx1() to DynamicOcpSize
				auto const nx_next = stage->B().rows();

				// Update D (which is initialized to 0)
				//K::middle_rows(Dc_, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs_.nu() - nu),
				//									  stage->C() * B, stage->D(), K::zero(sz.nc(), cs_.nu() - nu - sz.nu());
				submatrix(Dc_, nc          ,  0, sz.nx(),      nu) = B_;
				submatrix(Dc_, nc + sz.nx(),  0, sz.nc(),      nu) = stage->C() * B_;
				submatrix(Dc_, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->D();

				// Update lbd
				subvector(lbd_, nc          , sz.nx()) = stage->lbx() - g[k];
				subvector(lbd_, nc + sz.nx(), sz.nc()) = stage->lbd() - stage->C() * g[k];

				// Update ubd
				subvector(ubd_, nc          , sz.nx()) = stage->ubx() - g[k];
				subvector(ubd_, nc + sz.nx(), sz.nc()) = stage->ubd() - stage->C() * g[k];

				// Uplate lbu, ubu
				subvector(lbu_, nu, sz.nu()) = stage->lbu();
				subvector(ubu_, nu, sz.nu()) = stage->ubu();

				// Update B
				{
					DynamicMatrix<Real> B_next(nx_next, nu + sz.nu());
					submatrix(B_next, 0,  0, nx_next,      nu) = stage->A() * B_;
					submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->B();
					B_ = std::move(B_next);
				}

				// Update indices
				nu += sz.nu();
				nc += sz.nx() + sz.nc();
			}

			// Calculate r by backward substitution
			{
				r(N - 1) = stages[N - 1].r() + trans(stages[N - 1].S()) * g[N - 1];
				DynamicVector<Real> w = stages[N - 1].q() + stages[N - 1].Q() * g[N - 1];

				for (size_t k = N - 1; k --> 0 ;)
				{
					// Update r
					r(k) = stages[k].r() + trans(stages[k].S()) * g[k] + trans(stages[k].B()) * w;
					w = stages[k].q() + stages[k].Q() * g[k] + trans(stages[k].A()) * w;
				}
			}

			// Set the state bounds of the condensed stage<tmpc/ocp/DynamicOcpSize.hpp>
			qp.lbx(first_->lbx());
			qp.ubx(first_->ubx());

			// Assign the output values.				
			qp.A(F[N]);
			qp.B(B_);
			qp.b(g[N]);				
			qp.C(Cc_);
			qp.D(Dc_);
			qp.lbd(lbd_);
			qp.ubd(ubd_);
			qp.lbu(lbu_);
			qp.ubu(ubu_);
			qp.Q(Qc_);
			qp.R(Rc_);
			qp.S(Sc_);
			qp.q(qc_);
			qp.r(rc_);

			// TODO: implement soft constraints matrices recalculation.
			qp.Zl(sNaN<Real>());
			qp.Zu(sNaN<Real>());
			qp.zl(sNaN<Real>());
			qp.zu(sNaN<Real>());

			// TODO: what happens to the soft constraints index?
		}


		template <typename InIter>
		CondensingN2(InIter const& sz_first, InIter const& sz_last, size_t nx_next = 0)
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
		// Init G.
		,	G_(sumNx(sz_first, sz_last) + nx_next, sumNu(sz_first, sz_last), 0)
		,	blockG_(G_, 
			blockRowSizeG(sz_first, sz_last, nx_next),
			make_transform_iterator_range(sz_first, sz_last, [] (DynamicOcpSize const& s) { return s.nu(); })
			)
		,	size_(sz_first, sz_last)
		{
			size_t const N = std::distance(sz_first, sz_last);
			//g_.resize(N);
			
			// Init cumulative sizes.
			cumNu_.reserve(N);
			cumNc_.reserve(N);

			size_t nu = 0;
			size_t nc = 0;

			for (auto sz = size_.begin(); sz != size_.end(); ++sz)
			{
				cumNu_.push_back(nu);
				cumNc_.push_back(nc);
				nu += sz->nu();

				if (sz != size_.begin())
					nc += sz->nx();
					
				nc += sz->nc();
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
		DynamicOcpSize const cs_;

		DynamicMatrix<Real> Qc_;
		DynamicMatrix<Real> Rc_;
		DynamicMatrix<Real> Sc_;
		DynamicVector<Real> qc_;
		DynamicVector<Real> rc_;

		DynamicMatrix<Real> Cc_;
		DynamicMatrix<Real> Dc_;
		DynamicVector<Real> lbd_;
		DynamicVector<Real> ubd_;

		DynamicVector<Real> lbu_;
		DynamicVector<Real> ubu_;

		DynamicMatrix<Real> B_;
		DynamicMatrix<Real> G_;
		BlockMatrixView<DynamicMatrix<Real>> blockG_;

		/*
		class PerStageData
		{
		public:
			PerStageData(CondensingN2& c, size_t nx, size_t nu, size_t nc, size_t ns)
			:	

		private:
			CondensingN2& 

			DynamicOcpSize stageSize_;
			DynamicOcpSize cumSize_;

			DynamicVector<Real> g_;
			Subvector<Kernel, DynamicVector<Real>> r_;
		};
		*/

		//std::vector<DynamicVector<Real>> g_;
		std::vector<DynamicOcpSize> size_;
		std::vector<size_t> cumNu_;
		std::vector<size_t> cumNc_;


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


		decltype(auto) Cx(size_t i)
		{
			assert(i > 0);
			return submatrix(Cc_, cumNc_[i], 0, size_[i].nx(), columns(Cc_));
		}


		decltype(auto) Cc(size_t i)
		{
			return submatrix(Cc_, cumNc_[i] + size_[i].nx(), 0, size_[i].nc(), columns(Cc_));
		}


		template <typename InputIterator>
		static auto blockRowSizeG(InputIterator sz_first, InputIterator sz_last, size_t nx_next)
		{
			std::vector<size_t> sx;
			sx.reserve(std::distance(sz_first, sz_last) + nx_next);
			
			while (sz_first != sz_last)
				sx.push_back(sz_first++->nx());
			sx.push_back(nx_next);
			
			return sx;
		}
	};
}
