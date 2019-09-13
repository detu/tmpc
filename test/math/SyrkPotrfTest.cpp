#include <tmpc/math/SyrkPotrf.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
    TEST(MathTest, testSyrkPotrf)
    {
        size_t const m = 5, k = 4;  // <-- Sic!

        // Init Blaze matrices
        //
        blaze::DynamicMatrix<double, blaze::columnMajor> A(k, m), C(m, m), D, D1;
        randomize(A);
        makePositiveDefinite(C);

        // Do syrk-potrf with Blaze
        //
        llh(C + trans(A) * A, D);
        // std::cout << "D=\n" << D;

        // Do syrk-potrf with tmpc
        syrkPotrf(A, C, D1);

        // Print the result from syrkPotrf
        // std::cout << "D1=\n" << D1;
        TMPC_EXPECT_APPROX_EQ(D1, D, 1e-10, 1e-10);
    }
}
