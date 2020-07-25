#include <tmpc/json/JsonQp.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    namespace
    {
        template <typename MT, bool SO>
        auto resize(blaze::Matrix<MT, SO>&& m, std::pair<size_t, size_t> dims)
        {
            (~m).resize(dims.first, dims.second, true);
            return ~m;
        }
    }


    TEST(JsonQpTest, testFromJson)
    {
        auto const j = R"(
        {
            "branching": [2, 0, 0],
            "nodes": [
                {
                    "nx": 3, "nu": 2, "nc": 2,
                    "Q": [
                        [10, 1, 0],
                        [1, 15, 2],
                        [0, 2, 20]
                    ],
                    "R": [
                        [11, 2],
                        [2, 12]
                    ],
                    "S": [
                        [0.1, 0.2, 0.3],
                        [0.4, 0.5, 0.6]
                    ],
                    "q": [1, 2, 3],
                    "r": [4, 5],
                    "lx": [-100, -200, -300],
                    "ux": [100, 200, 300],
                    "lu": [-1.1, -2.2],
                    "uu": [1.1, 2.2],
                    "C": [
                        [2.1, 3.1, 4.1],
                        [5.1, 6.1, 7.1]
                    ],
                    "D": [
                        [2.2, 3.2],
                        [5.2, 6.2]
                    ],
                    "ld": [-20, -30],
                    "ud": [20, 30],
                    "zl": [],
                    "zu": [],
                    "Zl": [],
                    "Zu": []
                },
                {
                    "nx": 2, "nu": 0, "nc": 1,
                    "Q": [
                        [15, 2],
                        [2, 20]
                    ],
                    "R": [],
                    "S": [],
                    "q": [2, 3],
                    "r": [],
                    "lx": [-100, -200],
                    "ux": [100, 200],
                    "lu": [],
                    "uu": [],
                    "C": [
                        [2.1, 4.1]
                    ],
                    "D": [],
                    "ld": [-30],
                    "ud": [30],
                    "zl": [],
                    "zu": [],
                    "Zl": [],
                    "Zu": []
                },
                {
                    "nx": 2, "nu": 0, "nc": 1,
                    "Q": [
                        [15, 2],
                        [2, 20]
                    ],
                    "R": [],
                    "S": [],
                    "q": [2, 3],
                    "r": [],
                    "lx": [-100, -200],
                    "ux": [100, 200],
                    "lu": [],
                    "uu": [],
                    "C": [
                        [2.1, 4.1]
                    ],
                    "D": [],
                    "ld": [-30],
                    "ud": [30],
                    "zl": [],
                    "zu": [],
                    "Zl": [],
                    "Zu": []
                }
            ],
            "edges": [
                {
                    "from": 0,
                    "to": 1,
                    "A": [
                        [1, 2, 3],
                        [4, 5, 6]
                    ],
                    "B": [
                        [7, 8],
                        [9, 10]
                    ],
                    "b": [11, 12]
                },
                {
                    "from": 0,
                    "to": 2,
                    "A": [
                        [-1, -2, -3],
                        [-4, -5, -6]
                    ],
                    "B": [
                        [-7, -8],
                        [-9, -10]
                    ],
                    "b": [-11, -12]
                }
            ]
        }
        )"_json;


        using Vec = blaze::DynamicVector<double>;
        using Mat = blaze::DynamicMatrix<double>;

        DynamicOcpQp<double> qp = j;

        EXPECT_EQ(num_vertices(qp.graph()), j["nodes"].size());
        EXPECT_EQ(num_edges(qp.graph()), j["edges"].size());
        EXPECT_EQ(qp.size(), (DynamicOcpSize {
            OcpTree {2, 0, 0}, 
            {{3, 2, 2, 0}, {2, 0, 1, 0}, {2, 0, 1, 0}}
        }));

        for (auto v : vertices(qp.graph()))
        {
            TMPC_EXPECT_EQ(qp.Q(v), j["nodes"][v]["Q"].get<Mat>());
            TMPC_EXPECT_EQ(qp.q(v), j["nodes"][v]["q"].get<Vec>());
            TMPC_EXPECT_EQ(qp.lx(v), j["nodes"][v]["lx"].get<Vec>());
            TMPC_EXPECT_EQ(qp.ux(v), j["nodes"][v]["ux"].get<Vec>());
            TMPC_EXPECT_EQ(qp.C(v), resize(j["nodes"][v]["C"].get<Mat>(), {qp.size().nc(v), qp.size().nx(v)}));
            TMPC_EXPECT_EQ(qp.ld(v), j["nodes"][v]["ld"].get<Vec>());
            TMPC_EXPECT_EQ(qp.ud(v), j["nodes"][v]["ud"].get<Vec>());
        }

        for (auto v : qp.graph().branchVertices())
        {
            TMPC_EXPECT_EQ(qp.R(v), j["nodes"][v]["R"].get<Mat>());
            TMPC_EXPECT_EQ(qp.S(v), resize(j["nodes"][v]["S"].get<Mat>(), {qp.size().nu(v), qp.size().nx(v)}));
            TMPC_EXPECT_EQ(qp.r(v), j["nodes"][v]["r"].get<Vec>());
            TMPC_EXPECT_EQ(qp.lu(v), j["nodes"][v]["lu"].get<Vec>());
            TMPC_EXPECT_EQ(qp.uu(v), j["nodes"][v]["uu"].get<Vec>());
            TMPC_EXPECT_EQ(qp.D(v), resize(j["nodes"][v]["D"].get<Mat>(), {qp.size().nc(v), qp.size().nu(v)}));
        }
        
        for (auto e : edges(qp.graph()))
        {
            // Check values
            TMPC_EXPECT_EQ(qp.A(e), j["edges"][e]["A"].get<Mat>());
            TMPC_EXPECT_EQ(qp.B(e), j["edges"][e]["B"].get<Mat>());
            TMPC_EXPECT_EQ(qp.b(e), j["edges"][e]["b"].get<Vec>());
        }
    }


    /// @brief Check that a QP remains the same after converting to JSON and back.
    ///
    TEST(JsonQpTest, testToFromJson)
    {
        using Vec = blaze::DynamicVector<double>;
        using Mat = blaze::DynamicMatrix<double>;

        OcpTree const g {3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
        DynamicOcpSize const size {g, 3, 2, 1};

        DynamicOcpQp<double> qp0 {size};
        randomize(qp0);

        json j = qp0;
        DynamicOcpQp<double> qp1 = j;

        ASSERT_EQ(qp1.graph(), qp0.graph());

        for (auto v : vertices(qp1.graph()))
        {
            TMPC_EXPECT_EQ(qp1.Q(v), qp0.Q(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.q(v), qp0.q(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.lx(v), qp0.lx(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.ux(v), qp0.ux(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.C(v), qp0.C(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.ld(v), qp0.ld(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.ud(v), qp0.ud(v)) << "at v=" << v;
        }

        for (auto v : qp1.graph().branchVertices())
        {
            TMPC_EXPECT_EQ(qp1.R(v), qp0.R(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.S(v), qp0.S(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.r(v), qp0.r(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.lu(v), qp0.lu(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.uu(v), qp0.uu(v)) << "at v=" << v;
            TMPC_EXPECT_EQ(qp1.D(v), qp0.D(v)) << "at v=" << v;
        }
        
        for (auto e : edges(qp1.graph()))
        {
            // Check values
            TMPC_EXPECT_EQ(qp1.A(e), qp0.A(e)) << "at e=" << e;
            TMPC_EXPECT_EQ(qp1.B(e), qp0.B(e)) << "at e=" << e;
            TMPC_EXPECT_EQ(qp1.b(e), qp0.b(e)) << "at e=" << e;
        }
    }
}
