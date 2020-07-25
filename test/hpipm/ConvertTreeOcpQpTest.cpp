#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/Randomize.hpp>

#include <tmpc/hpipm/ConvertTreeOcpQp.hpp>

#include <tmpc/Testing.hpp>

#include <blasfeo/BlazeInterop.hpp>


namespace tmpc :: testing
{
	TEST(ConvertTreeOcpQpTest, testTree)
	{
		using Real = double;

		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;
		size_t constexpr NT = 2;

		OcpTree const g {
			3,
			2, 2, 2,
			1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0
		};
		// OcpTree const g {
		// 	2,
		// 	0, 0
		// };
		DynamicOcpSize const size {g, NX, NU, NC, 0, NCT, false};

		DynamicOcpQp<Real> qp {size};
		hpipm::TreeOcpQpDim<Real> dim {size};
		hpipm::TreeOcpQp<Real> hpipm_qp {dim};
		
		randomize(qp);
		convertQp(qp, hpipm_qp);

		blaze::DynamicMatrix<Real> M;
		blaze::DynamicVector<Real> v;

		auto const * const nu = hpipm_qp.dim->nu;
		auto const * const nx = hpipm_qp.dim->nx;
		auto const * const ng = hpipm_qp.dim->ng;

		for (auto ii : vertices(qp.graph()))
		{
			// Checking RSQrq:
			//
			// CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, c_ptr);
			//
			blasfeo::unpack(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], hpipm_qp.RSQrq[ii], 0, 0, M);

			if (out_degree(ii, qp.graph()) > 0)
			{
				// std::cout << "R=\n" << qp.R(ii);
				// std::cout << "S'=\n" << trans(qp.S(ii));
				// std::cout << "r'=\n" << trans(qp.r(ii));
				TMPC_EXPECT_EQ(submatrix(M, 0, 0, nu[ii], nu[ii]), qp.R(ii));
				TMPC_EXPECT_EQ(submatrix(M, nu[ii], 0, nx[ii], nu[ii]), trans(qp.S(ii)));
				TMPC_EXPECT_EQ(subvector(row(M, nu[ii] + nx[ii]), 0, nu[ii]), trans(qp.r(ii)));			
			}

			// std::cout << "Q=\n" << qp.Q(ii);
			// std::cout << "q'=\n" << trans(qp.q(ii));			
			TMPC_EXPECT_EQ(submatrix(M, nu[ii], nu[ii], nx[ii], nx[ii]), qp.Q(ii));
			TMPC_EXPECT_EQ(subvector(row(M, nu[ii] + nx[ii]), nu[ii], nx[ii]), trans(qp.q(ii)));

			// std::cout << "RSQrq=\n" << M;

			// Checking DCt:
			//
			// CREATE_STRMAT(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, c_ptr);
			//
			blasfeo::unpack(nu[ii] + nx[ii], ng[ii], hpipm_qp.DCt[ii], 0, 0, M);

			if (out_degree(ii, qp.graph()) > 0)
			{
				// std::cout << "D'=\n" << trans(qp.D(ii));
				TMPC_EXPECT_EQ(submatrix(M, 0, 0, nu[ii], ng[ii]), trans(qp.D(ii)));
			}

			// std::cout << "C'=\n" << trans(qp.C(ii));
			TMPC_EXPECT_EQ(submatrix(M, nu[ii], 0, nx[ii], ng[ii]), trans(qp.C(ii)));
		}


		// Checking BAbt:
		//
		// CREATE_STRMAT(nu[idxdad]+nx[idxdad]+1, nx[idx], qp->BAbt+ii, c_ptr);
		//
		for (auto e : edges(qp.graph()))
		{
			auto const idx = target(e, qp.graph());
			auto const idxdad = source(e, qp.graph());

			blasfeo::unpack(nu[idxdad] + nx[idxdad] + 1, nx[idx], hpipm_qp.BAbt[e], 0, 0, M);

			// std::cout << "B'=\n" << trans(qp.B(e));
			// std::cout << "A'=\n" << trans(qp.A(e));
			// std::cout << "b'=\n" << trans(qp.b(e));
			// std::cout << "BAbt=\n" << M;
			
			TMPC_EXPECT_EQ(submatrix(M, 0, 0, nu[idxdad], nx[idx]), trans(qp.B(e))) << " at edge " << e;
			TMPC_EXPECT_EQ(submatrix(M, nu[idxdad], 0, nx[idxdad], nx[idx]), trans(qp.A(e))) << " at edge " << e;
			TMPC_EXPECT_EQ(row(M, nu[idxdad] + nx[idxdad]), trans(qp.b(e))) << " at edge " << e;
		}
	}
}