#include <tmpc/json/JsonQp.hpp>
#include <tmpc/json/JsonBlaze.hpp>

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/core/Range.hpp>

#include <gtest/gtest.h>


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


    TEST(JsonQpTest, test_ocpQpToFromJson)
    {
        using Kernel = BlazeKernel<double>;
        using Stage = OcpQp<Kernel>;

        Stage stage0 {OcpSize {3, 2, 0}, 2};
        Stage stage1 {OcpSize {2, 1, 0}, 0};
        
        randomize(stage0);
        randomize(stage1);

        std::vector<Stage> qp;
        qp.push_back(stage0);
        qp.push_back(stage1);

        json j = ocpQpToJson(qp);

        std::cout << std::setw(4) << j << std::endl;
    }


    TEST(JsonQpTest, test_ctorFromJson)
    {
        auto const j = R"(
        {
            "nodes": [
                {
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
                        [0.1, 0.2],
                        [0.3, 0.4],
                        [0.5, 0.6]
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


        using K = BlazeKernel<double>;
        using Vec = DynamicVector<K>;
        using Mat = DynamicMatrix<K>;

        JsonQp<K> json_qp {j};

        EXPECT_EQ(num_vertices(json_qp.graph()), j["nodes"].size());
        EXPECT_EQ(num_edges(json_qp.graph()), j["edges"].size());
        EXPECT_EQ(get(json_qp.size(), 0), OcpSize(3, 2, 2, 0));
        EXPECT_EQ(get(json_qp.size(), 1), OcpSize(2, 0, 1, 0));

        for (auto v : make_iterator_range(vertices(json_qp.graph())))
        {
            EXPECT_EQ(get(json_qp.Q(), v), Mat(j["nodes"][v]["Q"]));
            EXPECT_EQ(get(json_qp.R(), v), Mat(j["nodes"][v]["R"]));
            EXPECT_EQ(get(json_qp.S(), v), resize(Mat(j["nodes"][v]["S"]), get(size_S(json_qp.size()), v)));
            EXPECT_EQ(get(json_qp.q(), v), Vec(j["nodes"][v]["q"]));
            EXPECT_EQ(get(json_qp.r(), v), Vec(j["nodes"][v]["r"]));

            EXPECT_EQ(get(json_qp.lx(), v), Vec(j["nodes"][v]["lx"]));
            EXPECT_EQ(get(json_qp.ux(), v), Vec(j["nodes"][v]["ux"]));
            EXPECT_EQ(get(json_qp.lu(), v), Vec(j["nodes"][v]["lu"]));
            EXPECT_EQ(get(json_qp.uu(), v), Vec(j["nodes"][v]["uu"]));

            EXPECT_EQ(get(json_qp.C(), v), resize(Mat(j["nodes"][v]["C"]), get(size_C(json_qp.size()), v)));
            EXPECT_EQ(get(json_qp.D(), v), resize(Mat(j["nodes"][v]["D"]), get(size_D(json_qp.size()), v)));
            EXPECT_EQ(get(json_qp.ld(), v), Vec(j["nodes"][v]["ld"]));
            EXPECT_EQ(get(json_qp.ud(), v), Vec(j["nodes"][v]["ud"]));
        }

        
        for (auto e : make_iterator_range(edges(json_qp.graph())))
        {
            // Find corresponding edge in json
            auto const j_edges = j["edges"];
            auto const u = source(e, json_qp.graph());
            auto const v = target(e, json_qp.graph());
            auto const j_edge = std::find_if(j_edges.begin(), j_edges.end(), [u, v] (auto je) 
            {
                return je["from"] == u && je["to"] == v;
            });

            ASSERT_NE(j_edge, j_edges.end());
            auto const edge_id = std::distance(j_edges.begin(), j_edge);

            // Check values
            EXPECT_EQ(get(json_qp.A(), e), Mat(j["edges"][edge_id]["A"]));
            EXPECT_EQ(get(json_qp.B(), e), Mat(j["edges"][edge_id]["B"]));
            EXPECT_EQ(get(json_qp.b(), e), Vec(j["edges"][edge_id]["b"]));
        }
    }
}
