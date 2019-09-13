#include <tmpc/blasfeo/Blasfeo.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>

#include <memory>


namespace tmpc :: testing
{
    static auto alignedAlloc(size_t bytes)
    {
        return std::unique_ptr<char[], decltype(&std::free)>(
            reinterpret_cast<char *>(std::aligned_alloc(0x40, bytes)), &std::free);
    }


    TEST(BlasfeoTest, testSyrkPotrf)
    {
        size_t const m = 5, k = 4;  // <-- Sic!

        // Init Blaze matrices
        //
        blaze::DynamicMatrix<double, blaze::columnMajor> blaze_A(m, k), blaze_C(m, m), blaze_D(m, m);
        randomize(blaze_A);
        makePositiveDefinite(blaze_C);

        // Do syrk-potrf with Blaze
        //
        llh(blaze_C + blaze_A * trans(blaze_A), blaze_D);
        // std::cout << "blaze_D=\n" << blaze_D;

        // Init BLASFEO matrices
        //
        blasfeo::DynamicMatrix<double> blasfeo_A(blaze_A), blasfeo_C(blaze_C), blasfeo_D(m, m);
        
        // Do syrk-potrf with BLASFEO
        //
        // void blasfeo_dsyrk_dpotrf_ln(int m, int k, 
        //  struct blasfeo_dmat *sA, int ai, int aj, 
        //  struct blasfeo_dmat *sB, int bi, int bj, 
        //  struct blasfeo_dmat *sC, int ci, int cj, 
        //  struct blasfeo_dmat *sD, int di, int dj)        
        blasfeo_dsyrk_dpotrf_ln(m, k, &blasfeo_A, 0, 0, &blasfeo_A, 0, 0, &blasfeo_C, 0, 0, &blasfeo_D, 0, 0);

        // Copy the resulting D matrix from BLASFEO to Blaze
        blaze::DynamicMatrix<double, blaze::columnMajor> blaze_blasfeo_D;
        blasfeo_D.unpack(blaze_blasfeo_D);

        // Print the result from BLASFEO
        // std::cout << "blaze_D=\n" << blaze_blasfeo_D;

        TMPC_EXPECT_APPROX_EQ(blaze_blasfeo_D, blaze_D, 1e-10, 1e-10);
    }
}
