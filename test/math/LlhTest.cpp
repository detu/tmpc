#include <tmpc/math/Llh.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    template <typename MatrixTypeA, typename MatrixTypeL>
    static void llhStaticTestImpl()
    {
        MatrixTypeA A;
        MatrixTypeL L, L_ref;
        makePositiveDefinite(A);

        blaze::llh(A, L_ref);
        tmpc::llh(A, L);

        TMPC_EXPECT_APPROX_EQ(L, L_ref, 1e-10, 0.);
    }


    template <typename MatrixTypeA>
    static void llhInplaceStaticTestImpl()
    {
        MatrixTypeA A;
        blaze::LowerMatrix<MatrixTypeA> L_ref;
        makePositiveDefinite(A);

        blaze::llh(A, L_ref);
        tmpc::llh(A);

        TMPC_EXPECT_APPROX_EQ(A, L_ref, 1e-10, 0.);
    }


    TEST(MathTest, test_Llh_double_1_1_StaticMatrix_rowMajor_StaticMatrix_rowMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 1, 1, blaze::rowMajor>, blaze::StaticMatrix<double, 1, 1, blaze::rowMajor>>();
    }


    TEST(MathTest, test_Llh_double_1_1_StaticMatrix_columnMajor_StaticMatrix_columnMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 1, 1, blaze::columnMajor>, blaze::StaticMatrix<double, 1, 1, blaze::columnMajor>>();
    }


    TEST(MathTest, test_Llh_double_2_2_StaticMatrix_rowMajor_StaticMatrix_rowMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 2, 2, blaze::rowMajor>, blaze::StaticMatrix<double, 2, 2, blaze::rowMajor>>();
    }


    TEST(MathTest, test_Llh_double_2_2_StaticMatrix_columnMajor_StaticMatrix_columnMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 2, 2, blaze::columnMajor>, blaze::StaticMatrix<double, 2, 2, blaze::columnMajor>>();
    }


    TEST(MathTest, test_Llh_double_10_10_StaticMatrix_rowMajor_StaticMatrix_rowMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 10, 10, blaze::rowMajor>, blaze::StaticMatrix<double, 10, 10, blaze::rowMajor>>();
    }


    TEST(MathTest, test_Llh_double_10_10_StaticMatrix_columnMajor_StaticMatrix_columnMajor)
    {
        llhStaticTestImpl<blaze::StaticMatrix<double, 10, 10, blaze::columnMajor>, blaze::StaticMatrix<double, 10, 10, blaze::columnMajor>>();
    }


    TEST(MathTest, test_Llh_double_10_10_SymmetricStaticMatrix_columnMajor_LowerStaticMatrix_columnMajor)
    {
        llhStaticTestImpl<
            blaze::SymmetricMatrix<blaze::StaticMatrix<double, 10, 10, blaze::columnMajor>>, 
            blaze::LowerMatrix<blaze::StaticMatrix<double, 10, 10, blaze::columnMajor>>
            >();
    }


    // ---------------------------------------------------------
    //
    // In-place LLH tests
    //
    // ---------------------------------------------------------

    TEST(MathTest, test_LlhInplace_double_1_1_StaticMatrix_rowMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 1, 1, blaze::rowMajor>>();
    }


    TEST(MathTest, test_LlhInplace_double_1_1_StaticMatrix_columnMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 1, 1, blaze::columnMajor>>();
    }


    TEST(MathTest, test_LlhInplace_double_2_2_StaticMatrix_rowMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 2, 2, blaze::rowMajor>>();
    }


    TEST(MathTest, test_LlhInplace_double_2_2_StaticMatrix_columnMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 2, 2, blaze::columnMajor>>();
    }


    TEST(MathTest, test_LlhInplace_double_10_10_StaticMatrix_rowMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 10, 10, blaze::rowMajor>>();
    }


    TEST(MathTest, test_LlhInplace_double_10_10_StaticMatrix_columnMajor)
    {
        llhInplaceStaticTestImpl<blaze::StaticMatrix<double, 10, 10, blaze::columnMajor>>();
    }


    // TEST(MathTest, testLlh_2x2_columnMajor)
    // {
    //     blaze::StaticMatrix<double, 2, 2, blaze::columnMajor> m {
    //         {1.1, 0.},
    //         {0., 2.2}
    //     };

    //     TMPC_EXPECT_APPROX_EQ(expm(m), (blaze::StaticMatrix<double, 2, 2> {
    //         {3.0042,         0.},
    //         {0.,    9.0250}
    //     }), 1e-4, 0.);
    // }


    // TEST(MathTest, testLlh_5x5_rowMajor)
    // {
    //     blaze::StaticMatrix<double, 5, 5, blaze::rowMajor> m {
    //         {0.814723686393179,   0.097540404999410,   0.157613081677548,   0.141886338627215,   0.655740699156587},
    //         {0.905791937075619,   0.278498218867048,   0.970592781760616,   0.421761282626275,   0.035711678574190},
    //         {0.126986816293506,   0.546881519204984,   0.957166948242946,   0.915735525189067,   0.849129305868777},
    //         {0.913375856139019,   0.957506835434298,   0.485375648722841,   0.792207329559554,   0.933993247757551},
    //         {0.632359246225410,   0.964888535199277,   0.800280468888800,   0.959492426392903,   0.678735154857773}
    //     };

    //     TMPC_EXPECT_APPROX_EQ(expm(m), (blaze::StaticMatrix<double, 5, 5> {
    //         {4.204256643154649,   2.029099147640609,   2.316816191101029,   2.301099526225243,   2.922876596587000},
    //         {4.053307051154341,   3.927899946799651,   4.139010062772191,   3.523504905263845,   3.285583505067035},
    //         {4.980810366851268,   5.016681688555301,   7.072138611807674,   5.965554741935332,   5.677295953715725},
    //         {6.384560022826461,   5.444182034921380,   5.724925285245762,   6.848391976712574,   5.966624393908855},
    //         {6.093780099340079,   5.591149359420235,   6.233313648084653,   6.208172974358048,   6.841880083058702}
    //     }), 1e-14, 0.);
    // }


    // TEST(MathTest, testLlh_5x5_columnMajor)
    // {
    //     blaze::StaticMatrix<double, 5, 5, blaze::columnMajor> m {
    //         {0.814723686393179,   0.097540404999410,   0.157613081677548,   0.141886338627215,   0.655740699156587},
    //         {0.905791937075619,   0.278498218867048,   0.970592781760616,   0.421761282626275,   0.035711678574190},
    //         {0.126986816293506,   0.546881519204984,   0.957166948242946,   0.915735525189067,   0.849129305868777},
    //         {0.913375856139019,   0.957506835434298,   0.485375648722841,   0.792207329559554,   0.933993247757551},
    //         {0.632359246225410,   0.964888535199277,   0.800280468888800,   0.959492426392903,   0.678735154857773}
    //     };

    //     TMPC_EXPECT_APPROX_EQ(expm(m), (blaze::StaticMatrix<double, 5, 5> {
    //         {4.204256643154649,   2.029099147640609,   2.316816191101029,   2.301099526225243,   2.922876596587000},
    //         {4.053307051154341,   3.927899946799651,   4.139010062772191,   3.523504905263845,   3.285583505067035},
    //         {4.980810366851268,   5.016681688555301,   7.072138611807674,   5.965554741935332,   5.677295953715725},
    //         {6.384560022826461,   5.444182034921380,   5.724925285245762,   6.848391976712574,   5.966624393908855},
    //         {6.093780099340079,   5.591149359420235,   6.233313648084653,   6.208172974358048,   6.841880083058702}
    //     }), 1e-14, 0.);
    // }
}