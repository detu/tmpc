/*
 * Unit tests for common functionality of QP workspaces.
 */

#include "gtest_tools_eigen.hpp"

#include <tmpc/qp/Printing.hpp>
#include <tmpc/qp/QpOasesWorkspace.hpp>

#ifdef WITH_HPMPC
#include <tmpc/qp/HpmpcWorkspace.hpp>
#endif

#ifdef WITH_HPIPM
#include <tmpc/qp/HpipmWorkspace.hpp>
#endif

#include <gtest/gtest.h>

#include <iostream>
#include <array>

using namespace tmpc;

template <typename WS>
class QpWorkspaceTest : public ::testing::Test
{
public:
	using Workspace = WS;

protected:
	using Real = typename Workspace::Real;
	using Vector = DynamicVector<Real>;
	using Matrix = DynamicMatrix<Real>;

	QpWorkspaceTest()
	:	size_{
			QpSize(2, 3, 4),
			QpSize(5, 6, 7),
			QpSize(8, 9, 10)
		}
	,	ws_(size_.begin(), size_.end())
	{
	}

	std::array<QpSize, 3> const size_;
	Workspace ws_;
};

using WsTypes = ::testing::Types<
	QpOasesWorkspace
	#ifdef WITH_HPMPC
	,	HpmpcWorkspace<double>
	#endif
	#ifdef WITH_HPIPM
	,	HpipmWorkspace<double>
	#endif
	>;

TYPED_TEST_CASE(QpWorkspaceTest, WsTypes);

TYPED_TEST(QpWorkspaceTest, testQpInterface)
{
	auto const N = this->size_.size();

	std::vector<typename TestFixture::Matrix> Q(N);
	std::vector<typename TestFixture::Vector> q(N);
	std::vector<typename TestFixture::Vector> r(N);

	std::vector<typename TestFixture::Matrix> S(N);
	std::vector<typename TestFixture::Matrix> R(N);
	std::vector<typename TestFixture::Matrix> A(N);
	std::vector<typename TestFixture::Matrix> B(N);
	std::vector<typename TestFixture::Vector> b(N);

	std::vector<typename TestFixture::Vector> x_min(N), x_max(N);
	std::vector<typename TestFixture::Vector> u_min(N), u_max(N);
	std::vector<typename TestFixture::Vector> d_min(N), d_max(N);

	// Writing random data
	Rand<typename TestFixture::Matrix> rand_matrix;
	Rand<typename TestFixture::Vector> rand_vector;

	for (std::size_t i = 0; i < N; ++i)
	{
		auto const& sz = this->size_[i];
		auto const nx1 = i + 1 < N ? this->size_[i + 1].nx() : 0;
		auto& stage = this->ws_[i];

		stage.Q(Q[i] = rand_matrix.generate(sz.nx(), sz.nx()));
		stage.R(R[i] = rand_matrix.generate(sz.nu(), sz.nu()));
		stage.S(S[i] = rand_matrix.generate(sz.nx(), sz.nu()));
		stage.q(q[i] = rand_vector.generate(sz.nx()));
		stage.r(r[i] = rand_vector.generate(sz.nu()));

		stage.A(A[i] = rand_matrix.generate(nx1, sz.nx()));
		stage.B(B[i] = rand_matrix.generate(nx1, sz.nu()));
		stage.b(b[i] = rand_vector.generate(nx1));

		stage.lbx(x_min[i] = rand_vector.generate(sz.nx()));
		stage.ubx(x_max[i] = rand_vector.generate(sz.nx()));
		stage.lbu(u_min[i] = rand_vector.generate(sz.nu()));
		stage.ubu(u_max[i] = rand_vector.generate(sz.nu()));
		stage.lbd(d_min[i] = rand_vector.generate(sz.nc()));
		stage.ubd(d_max[i] = rand_vector.generate(sz.nc()));
	}

	// Reading the data and checking that they are the same that we wrote
	for (std::size_t i = 0; i < N; ++i)
	{
		auto const& stage = this->ws_[i];

		EXPECT_EQ(print_wrap(stage.Q()), print_wrap(Q[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(stage.R()), print_wrap(R[i])) << "at i=" << i;
		EXPECT_EQ(print_wrap(stage.S()), print_wrap(S[i]));
		EXPECT_EQ(print_wrap(stage.q()), print_wrap(q[i]));
		EXPECT_EQ(print_wrap(stage.r()), print_wrap(r[i]));

		EXPECT_EQ(print_wrap(stage.A()), print_wrap(A[i]));
		EXPECT_EQ(print_wrap(stage.B()), print_wrap(B[i]));
		EXPECT_EQ(print_wrap(stage.b()), print_wrap(b[i]));

		EXPECT_EQ(print_wrap(stage.lbx()), print_wrap(x_min[i]));
		EXPECT_EQ(print_wrap(stage.ubx()), print_wrap(x_max[i]));
		EXPECT_EQ(print_wrap(stage.lbu()), print_wrap(u_min[i]));
		EXPECT_EQ(print_wrap(stage.ubu()), print_wrap(u_max[i]));
		EXPECT_EQ(print_wrap(stage.lbd()), print_wrap(d_min[i]));
		EXPECT_EQ(print_wrap(stage.ubd()), print_wrap(d_max[i]));
	}
}

TYPED_TEST(QpWorkspaceTest, testMatrixSizesCorrect)
{
	// Define dimensions
	unsigned constexpr NX = 2;
	unsigned constexpr NU = 1;
	unsigned constexpr NC = 0;
	unsigned constexpr NCT = 0;
	unsigned constexpr NT = 2;

	std::vector<QpSize> sz;
	sz.reserve(NT + 1);
	
	for (size_t i = 0; i < NT; ++i)
		sz.emplace_back(NX, NU, NC);
	sz.emplace_back(NX, 0, NCT);
	
	typename TestFixture::Workspace ws(sz.begin(), sz.end());

	for (std::size_t i = 0; i < sz.size(); ++i)
	{
		auto const& s = sz[i];
		auto const nx1 = i + 1 < sz.size() ? sz[i + 1].nx() : 0;
		auto const& stage = ws[i];

		EXPECT_EQ(rows   (stage.Q()), s.nx());
		EXPECT_EQ(columns(stage.Q()), s.nx());
		EXPECT_EQ(rows   (stage.R()), s.nu());
		EXPECT_EQ(columns(stage.R()), s.nu());
		EXPECT_EQ(rows   (stage.S()), s.nx());
		EXPECT_EQ(columns(stage.S()), s.nu());
		EXPECT_EQ(size   (stage.q()), s.nx());
		EXPECT_EQ(size   (stage.r()), s.nu());

		EXPECT_EQ(rows   (stage.A()),   nx1 );
		EXPECT_EQ(columns(stage.A()), s.nx());
		EXPECT_EQ(rows   (stage.B()),   nx1 );
		EXPECT_EQ(columns(stage.B()), s.nu());
		EXPECT_EQ(size   (stage.b()),   nx1 );

		EXPECT_EQ(size(stage.lbx()), s.nx());
		EXPECT_EQ(size(stage.ubx()), s.nx());
		EXPECT_EQ(size(stage.lbu()), s.nu());
		EXPECT_EQ(size(stage.ubu()), s.nu());
		EXPECT_EQ(size(stage.lbd()), s.nc());
		EXPECT_EQ(size(stage.ubd()), s.nc());
	}
}
