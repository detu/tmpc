/*
@brief Support convertion between OCP QP and JSON.
*/

#pragma once

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/json/Json.hpp>
#include <tmpc/json/JsonBlaze.hpp>
#include <tmpc/Exception.hpp>

#include <string>
#include <ranges>


namespace nlohmann 
{
    template <tmpc::OcpQp Qp>
    struct adl_serializer<Qp> 
    {
        // Here's the catch! You must provide a to_json method! Otherwise you
        // will not be able to convert move_only_type to json, since you fully
        // specialized adl_serializer on that type
        static void to_json(json& j, Qp const& qp) 
        {
            json j_nodes;
            json j_branching;

            for (auto v : vertices(qp.graph()))
            {
                json j_node {
                    {"nx", qp.size().nx(v)},
                    {"nu", qp.size().nu(v)},
                    {"nc", qp.size().nc(v)},
                    {"Q", qp.Q(v)},
                    {"q", qp.q(v)},
                    {"C", qp.C(v)},
                    {"ld", qp.ld(v)},
                    {"ud", qp.ud(v)},
                    {"lx", qp.lx(v)},
                    {"ux", qp.ux(v)},
                    // {"Zl", get(qp.Zl(), v)},
                    // {"Zu", get(qp.Zu(), v)},
                    // {"zl", get(qp.zl(), v)},
                    // {"zu", get(qp.zu(), v)}
                };

                if (out_degree(v, qp.graph()) > 0)
                {
                    j_node["R"] = qp.R(v);
                    j_node["S"] = qp.S(v);
                    j_node["r"] = qp.r(v);
                    j_node["D"] = qp.D(v);
                    j_node["lu"] = qp.lu(v);
                    j_node["uu"] = qp.uu(v);
                }

                j_nodes.push_back(j_node);
                j_branching.push_back(out_degree(v, qp.graph()));
            }

            json j_edges;

            for (auto e : edges(qp.graph()))
            {
                json j_edge {
                    {"A", qp.A(e)},
                    {"B", qp.B(e)},
                    {"b", qp.b(e)},
                    {"from", source(e, qp.graph())},
                    {"to", target(e, qp.graph())}
                };

                j_edges.push_back(j_edge);
            }

            j["branching"] = j_branching;
            j["nodes"] = j_nodes;
            j["edges"] = j_edges;
        }


        static Qp from_json(const json& j) 
        {
            tmpc::OcpTree g {j.at("branching") 
                | std::views::transform([] (auto const& n) -> std::size_t { return n; })};

            auto const j_nodes = j.at("nodes");
            auto const j_edges = j.at("edges");

            tmpc::DynamicOcpSize size {g, 
                [&j_nodes] (tmpc::OcpVertex v, tmpc::OcpTree const& g)
                {
                    return tmpc::OcpVertexSize {
                        j_nodes.at(v).at("nx"),
                        j_nodes.at(v).at("nu"),
                        j_nodes.at(v).at("nc")};
                }
            };

            Qp qp {size};

            for (auto v : vertices(g))
            {
                auto const j_node = j_nodes.at(v);
                
                qp.Q(v, reshapeIfNeeded(j_node.at("Q").get<Matrix>(), size.nx(v), size.nx(v)));
                qp.C(v, reshapeIfNeeded(j_node.at("C").get<Matrix>(), size.nc(v), size.nx(v)));
                qp.q(v, j_node.at("q").get<Vector>());
                qp.lx(v, j_node.at("lx").get<Vector>());
                qp.ux(v, j_node.at("ux").get<Vector>());
                qp.ld(v, j_node.at("ld").get<Vector>());
                qp.ud(v, j_node.at("ud").get<Vector>());
            }

            for (auto v : g.branchVertices())
            {
                auto const j_node = j_nodes.at(v);
                
                qp.R(v, reshapeIfNeeded(j_node.at("R").get<Matrix>(), size.nu(v), size.nu(v)));
                qp.S(v, reshapeIfNeeded(j_node.at("S").get<Matrix>(), size.nu(v), size.nx(v)));
                qp.D(v, reshapeIfNeeded(j_node.at("D").get<Matrix>(), size.nc(v), size.nu(v)));
                qp.r(v, j_node.at("r").get<Vector>());
                qp.lu(v, j_node.at("lu").get<Vector>());
                qp.uu(v, j_node.at("uu").get<Vector>());
            }

            for (auto e : edges(g))
            {
                auto const j_edge = j_edges.at(e);
                auto const u = source(e, qp.graph());
                auto const v = target(e, qp.graph());

                qp.A(e, reshapeIfNeeded(j_edge.at("A").get<Matrix>(), size.nx(v), size.nx(u)));
                qp.B(e, reshapeIfNeeded(j_edge.at("B").get<Matrix>(), size.nx(v), size.nu(u)));
                qp.b(e, j_edge.at("b").get<Vector>());
            }

            return qp;
        }


    private:
        using Real = typename Qp::Real;
        using Vector = blaze::DynamicVector<Real>;
        using Matrix = blaze::DynamicMatrix<Real>;

            
        static Matrix reshapeIfNeeded(Matrix&& val, size_t expected_rows, size_t expected_columns)
        {
            // If both val and expected_dim are empty, resize val to match expected_dim.
            if (isEmpty(val) && (expected_rows == 0 || expected_columns == 0))
                val.resize(expected_rows, expected_columns);

            // If val is Nx1 and expected_dim is 1xN, transpose val to match expected_dim.
            if (columns(val) == 1 && expected_rows == 1 && expected_columns == rows(val))
                transpose(val);

            return val;
        }
    };
}
