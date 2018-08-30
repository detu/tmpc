/*
@brief Support convertion between OCP QP and JSON.
*/

#pragma once

#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/core/Range.hpp>

#include "Json.hpp"

#include <string>


namespace tmpc
{
    template <typename QP>
    void to_json(json& j, OcpQpBase<QP> const& qp) 
    {
        j = json {
            {"Q", qp.Q()},
            {"R", qp.R()},
            {"S", qp.S()},
            {"q", qp.q()},
            {"r", qp.r()},
            {"C", qp.C()},
            {"D", qp.D()},
            {"ld", qp.lbd()},
            {"ud", qp.ubd()},
            {"Zl", qp.Zl()},
            {"Zu", qp.Zu()},
            {"zl", qp.zl()},
            {"zu", qp.zu()},
            {"idxs", qp.idxs()}
        };
    }


    template <typename QP>
    void from_json(json const& j, OcpQpVertexBase<QP>& qp)
    {
        using Kernel = typename QP::Kernel;
        using Vector = DynamicVector<Kernel>;
        using Matrix = DynamicMatrix<Kernel>;

        qp.Q(j["Q"].get<Matrix>());
        qp.R(j["R"].get<Matrix>());
        qp.S(j["S"].get<Matrix>());
        qp.q(j["q"].get<Vector>());
        qp.r(j["r"].get<Vector>());
        qp.Zl(j["Zl"].get<Matrix>());
        qp.Zu(j["Zu"].get<Matrix>());
        qp.zl(j["zl"].get<Vector>());
        qp.zu(j["zu"].get<Vector>());
        qp.lbx(j["lx"].get<Vector>());
        qp.ubx(j["ux"].get<Vector>());
        qp.lbu(j["lu"].get<Vector>());
        qp.ubu(j["uu"].get<Vector>());
        qp.idxs(j["idxs"]);
    }


    template <typename QP>
    void from_json(json const& j, OcpQpEdgeBase<QP>& qp)
    {
        using Kernel = typename QP::Kernel;
        using Vector = DynamicVector<Kernel>;
        using Matrix = DynamicMatrix<Kernel>;

        qp.A(j["A"].get<Matrix>());
        qp.B(j["B"].get<Matrix>());
        qp.b(j["b"].get<Vector>());
    }


    template <typename IteratorRange>
    json ocpQpToJson(IteratorRange const& qp)
    {
        json nodes = json::array();
        json edges = json::array();
        
        size_t node_index = 0;

        for (auto stage = qp.begin(); stage != qp.end(); ++stage, ++node_index)
        {
            nodes.push_back(json {
                {"Q", stage->Q()},
                {"R", stage->R()},
                {"S", stage->S()},
                {"q", stage->q()},
                {"r", stage->r()},
                {"C", stage->C()},
                {"D", stage->D()},
                {"ld", stage->lbd()},
                {"ud", stage->ubd()},
                {"Zl", stage->Zl()},
                {"Zu", stage->Zu()},
                {"zl", stage->zl()},
                {"zu", stage->zu()},
                {"lx", stage->lbx()},
                {"ux", stage->ubx()},
                {"lu", stage->lbu()},
                {"uu", stage->ubu()},
                {"idxs", stage->idxs()}
            });

            if (stage + 1 != qp.end())
            {
                edges.push_back(json {
                    {"from", node_index},
                    {"to", node_index + 1},
                    {"A", stage->A()},
                    {"B", stage->B()},
                    {"b", stage->b()}
                });
            }
        };

        return json {
            {"nodes", nodes},
            {"edges", edges}
        };
    }


    namespace detail
    {
        template <typename Descriptor>
        struct DescriptorTraits;

        
        template <>
        struct DescriptorTraits<OcpVertexDescriptor>
        {
            static char const * dictionaryKey()
            {
                return "nodes";
            }

            using index_property_t = boost::vertex_index_t;
        };


        template <>
        struct DescriptorTraits<OcpEdgeDescriptor>
        {
            static char const * dictionaryKey()
            {
                return "edges";
            }

            using index_property_t = boost::edge_index_t;
        };


        template <
            typename Key, 
            typename Value,
            typename IndexMap
        >
        class JsonQpPropertyMap
        {
        public:
            using key_type = Key; 
            using value_type = Value;
            using category = read_write_property_map_tag;

            JsonQpPropertyMap(json& j, std::string const& name, IndexMap index_map) 
            :   json_ {j}
            ,   name_ {name}
            ,   indexMap_ {index_map}
            {
            }


            friend void put(JsonQpPropertyMap const& m, Key const& key, Value const& val)
            {
                m.json_.at(get(m.indexMap_, key))[m.name_] = val;
            }


            friend Value get(JsonQpPropertyMap const& m, Key const& key)
            {
                auto const j_obj = m.json_.at(get(m.indexMap_, key));
                auto const j_val = j_obj.find(m.name_);

                // Return empty value if the m.name_ key is not present in json object.
                return j_val != j_obj.end() ? Value(j_val.value()) : Value();
            }


        private:
            json& json_;
            std::string const name_;
            IndexMap indexMap_;
        };



        template <
            typename Key, 
            typename Value,
            typename IndexMap
        >
        inline auto makeJsonQpPropertyMap(json& j, std::string const& name, IndexMap index_map)
        {
            return JsonQpPropertyMap<Key, Value, IndexMap>(j, name, index_map);
        }
    }


    template <typename Kernel>
    class JsonQp
    {
    public:
        using vertex_descriptor = OcpVertexDescriptor;
        using edge_descriptor = OcpEdgeDescriptor;


        JsonQp(tmpc::json j)
        :   json_(j)
        {
            auto const j_nodes = j.at("nodes");

            // Build OCP graph from json.
            size_t const N = j_nodes.size();
            graph_ = OcpGraph(N);

            // Add edges, with edge_index attribute equal to the index of an edge in "edges" json array.
            for (auto j_edge : j.at("edges") | indexed())
                add_edge(j_edge.value().at("from"), j_edge.value().at("to"), j_edge.index(), graph_);

            // Fill node sizes.
            size_.resize(N);
            std::transform(j_nodes.begin(), j_nodes.end(), size_.begin(), [] (auto const& j_v)
            {
                return OcpSize {
                    j_v.count("q") ? j_v["q"].size() : 0, 
                    j_v.count("r") ? j_v["r"].size() : 0, 
                    j_v.count("ld") ? j_v["ld"].size() : 0, 
                    j_v.count("zl") ? j_v["zl"].size() : 0};
            });
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto const& get_json() const
        {
            return json_;
        }


        auto size() const
        {
            return iterator_property_map(size_.begin(), get(vertex_index, graph_));
        }


        auto Q()
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "Q", get(vertex_index, graph_));
        }


        auto R()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_["nodes"], "R", get(vertex_index, graph_));
        }


        auto S()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_["nodes"], "S", get(vertex_index, graph_));
        }


        auto q()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "q", get(vertex_index, graph_));
        }


        auto r()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "r", get(vertex_index, graph_));
        }


        auto lx()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "lx", get(vertex_index, graph_));
        }


        auto ux()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "ux", get(vertex_index, graph_));
        }


        auto lu()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "lu", get(vertex_index, graph_));
        }


        auto uu()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "uu", get(vertex_index, graph_));
        }


        auto C()
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "C", get(vertex_index, graph_));
        }


        auto D()
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "D", get(vertex_index, graph_));
        }


        auto ld()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "ld", get(vertex_index, graph_));
        }


        auto ud()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_["nodes"], "ud", get(vertex_index, graph_));
        }


        auto A()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "A", get(edge_index, graph_));
        }


        auto B()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "B", get(edge_index, graph_));
        }


        auto b()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicVector<Kernel>>(
                json_["edges"], "b", get(edge_index, graph_));
        }


    private:
        // OCP graph
        OcpGraph graph_;

        // OCP node sizes
        std::vector<OcpSize> size_;

        // OCP JSON
        tmpc::json json_;
    };
}
