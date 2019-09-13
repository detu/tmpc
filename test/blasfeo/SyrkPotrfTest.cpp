#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

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
        auto mem_A = alignedAlloc(blasfeo_memsize_dmat(m, k));
        auto mem_C = alignedAlloc(blasfeo_memsize_dmat(m, m));
        auto mem_D = alignedAlloc(blasfeo_memsize_dmat(m, m));

        blasfeo_dmat blasfeo_A, blasfeo_C, blasfeo_D;
        blasfeo_create_dmat(m, k, &blasfeo_A, mem_A.get());
        blasfeo_create_dmat(m, m, &blasfeo_C, mem_C.get());
        blasfeo_create_dmat(m, m, &blasfeo_D, mem_D.get());

        // Copy Blaze matrices to BLASFEO matrices
        //
        // blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);
        blasfeo_pack_dmat(m, k, data(blaze_A), spacing(blaze_A), &blasfeo_A, 0, 0);
        blasfeo_pack_dmat(m, m, data(blaze_C), spacing(blaze_C), &blasfeo_C, 0, 0);
        
        // Do syrk-potrf with BLASFEO
        //
        // void blasfeo_dsyrk_dpotrf_ln(int m, int k, 
        //  struct blasfeo_dmat *sA, int ai, int aj, 
        //  struct blasfeo_dmat *sB, int bi, int bj, 
        //  struct blasfeo_dmat *sC, int ci, int cj, 
        //  struct blasfeo_dmat *sD, int di, int dj)        
        blasfeo_dsyrk_dpotrf_ln(m, k, &blasfeo_A, 0, 0, &blasfeo_A, 0, 0, &blasfeo_C, 0, 0, &blasfeo_D, 0, 0);

        // Copy the resulting D matrix from BLASFEO to Blaze
        //
        // void blasfeo_unpack_dmat(int m, int n, struct blasfeo_dmat *sA, int ai, int aj, double *B, int ldb);
        blaze::DynamicMatrix<double, blaze::columnMajor> blaze_blasfeo_D(m, m);
        blasfeo_unpack_dmat(m, m, &blasfeo_D, 0, 0, data(blaze_blasfeo_D), spacing(blaze_blasfeo_D));

        // Print the result from BLASFEO
        // std::cout << "blaze_D=\n" << blaze_blasfeo_D;

        TMPC_EXPECT_APPROX_EQ(blaze_blasfeo_D, blaze_D, 1e-10, 1e-10);
    }
}
