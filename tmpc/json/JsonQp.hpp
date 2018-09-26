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
#include <tmpc/Traits.hpp>

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
        template <typename Matrix>
        inline void reshapeIfNeeded(Matrix& val, std::pair<size_t, size_t> expected_dim)
        {
            // If both val and expected_dim are empty, resize val to match expected_dim.
            if (isEmpty(val) && (expected_dim.first == 0 || expected_dim.second == 0))
                val.resize(expected_dim.first, expected_dim.second);

            // If val is Nx1 and expected_dim is 1xN, transpose val to match expected_dim.
            if (columns(val) == 1 && expected_dim.first == 1 && expected_dim.second == rows(val))
                transpose(val);
        }


        template <typename Vector>
        inline void reshapeIfNeeded(Vector& val, size_t expected_dim)
        {
            // do nothing
        }


        template <typename Vector>
        inline void resize(Vector& val, size_t n)
        {
            val.resize(n);
        }


        template <typename Matrix>
        inline void resize(Matrix& val, std::pair<size_t, size_t> sz)
        {
            val.resize(sz.first, sz.second);
        }


        template <
            typename Key, 
            typename Value,
            typename IndexMap,
            typename SizeMap,
            bool Writeable = true
        >
        class JsonQpPropertyMap
        {
        public:
            using key_type = Key; 
            using value_type = Value;
            using category = std::conditional_t<Writeable, read_write_property_map_tag, readable_property_map_tag>;
            using JsonRef = std::conditional_t<Writeable, json&, json const&>;
            

            JsonQpPropertyMap(JsonRef j, std::string const& name, IndexMap index_map, SizeMap size_map) 
            :   json_ {j}
            ,   name_ {name}
            ,   indexMap_ {index_map}
            ,   sizeMap_ {size_map}
            {
            }


            friend std::enable_if_t<Writeable, void> put(JsonQpPropertyMap const& m, Key const& key, Value const& val)
            {
                if (dimensions(val) != get(m.sizeMap_, key))
                    throw std::invalid_argument("Invalid matrix/vector dimensions detected while trying to set JsonQp property \"" + m.name_ + "\"");

                m.json_.at(get(m.indexMap_, key))[m.name_] = val;
            }


            friend Value get(JsonQpPropertyMap const& m, Key const& key)
            {
                auto const j_obj = m.json_.at(get(m.indexMap_, key));
                auto const j_val = j_obj.find(m.name_);
                auto const expected_dim = get(m.sizeMap_, key);

                Value val;

                if (j_val != j_obj.end())
                {
                    from_json(j_val.value(), val);
                }
                else
                {
                    // If one of lx, lu, ux, uu is not present, return +-inf of appropriate size.
                    if (m.name_ == "lx" || m.name_ == "lu")
                    {
                        resize(val, expected_dim);
                        val = -inf<Scalar>();
                    }
                    else if (m.name_ == "ux" || m.name_ == "uu")
                    {
                        resize(val, expected_dim);
                        val = inf<Scalar>();
                    }

                    // Return empty value if the m.name_ key is not present in json object.
                }

                reshapeIfNeeded(val, expected_dim);

                if (dimensions(val) != expected_dim)
                    throw std::invalid_argument("Invalid matrix/vector dimensions detected while trying to get JsonQp property \"" + m.name_ + "\"");

                return val;
            }


            template <typename Range>
            void defaultInit(Range const& range)
            {
                std::vector<value_type> values(size(range));

                for (auto key : range)
                {
                    value_type val;
                    resize(val, get(sizeMap_, key));
                    val = 0.;

                    values.at(get(indexMap_, key)) = val;
                    //j_tmp.at(get(indexMap_, key))[name_] = val;
                }

                for (auto v : boost::adaptors::index(values))
                    json_[v.index()][name_] = v.value();
            }


        private:
            using Scalar = typename Value::ElementType;

            JsonRef json_;
            std::string const name_;
            IndexMap indexMap_;
            SizeMap sizeMap_;
        };



        template <
            typename Key, 
            typename Value,
            typename IndexMap,
            typename SizeMap
        >
        inline auto makeJsonQpPropertyMap(json& j, std::string const& name, IndexMap index_map, SizeMap size_map)
        {
            return JsonQpPropertyMap<Key, Value, IndexMap, SizeMap, true>(j, name, index_map, size_map);
        }


        template <
            typename Key, 
            typename Value,
            typename IndexMap,
            typename SizeMap
        >
        inline auto makeJsonQpPropertyMap(json const& j, std::string const& name, IndexMap index_map, SizeMap size_map)
        {
            return JsonQpPropertyMap<Key, Value, IndexMap, SizeMap, false>(j, name, index_map, size_map);
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


        template <typename SizeMap>
        JsonQp(OcpGraph const& g, SizeMap size_map)
        :   graph_{g}
        ,   size_(num_vertices(g))
        ,   json_{
            {"nodes", json::array()}, 
            {"edges", json::array()}
            }
        {
            auto const vert = vertices(g);
            copyProperty(size_map, iterator_property_map(size_.begin(), vertexIndex(g)), vert);

            Q().defaultInit(vert);
            R().defaultInit(vert);
            S().defaultInit(vert);
            q().defaultInit(vert);
            r().defaultInit(vert);
            lx().defaultInit(vert);
            ux().defaultInit(vert);
            lu().defaultInit(vert);
            uu().defaultInit(vert);
            C().defaultInit(vert);
            D().defaultInit(vert);
            ld().defaultInit(vert);
            ud().defaultInit(vert);

            auto const edg = edges(g);
            A().defaultInit(edg);
            B().defaultInit(edg);
            b().defaultInit(edg);
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
                json_.at("nodes"), "Q", get(vertex_index, graph_), size_Q(size()));
        }


        auto Q() const
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "Q", get(vertex_index, graph_), size_Q(size()));
        }


        auto R()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "R", get(vertex_index, graph_), size_R(size()));
        }


        auto R() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "R", get(vertex_index, graph_), size_R(size()));
        }


        auto S()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "S", get(vertex_index, graph_), size_S(size()));
        }


        auto S() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "S", get(vertex_index, graph_), size_S(size()));
        }


        auto q()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "q", get(vertex_index, graph_), size_x(size()));
        }


        auto q() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "q", get(vertex_index, graph_), size_x(size()));
        }


        auto r()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "r", get(vertex_index, graph_), size_u(size()));
        }


        auto r() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "r", get(vertex_index, graph_), size_u(size()));
        }


        auto lx()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "lx", get(vertex_index, graph_), size_x(size()));
        }


        auto lx() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "lx", get(vertex_index, graph_), size_x(size()));
        }


        auto ux()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ux", get(vertex_index, graph_), size_x(size()));
        }


        auto ux() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ux", get(vertex_index, graph_), size_x(size()));
        }


        auto lu()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "lu", get(vertex_index, graph_), size_u(size()));
        }


        auto lu() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "lu", get(vertex_index, graph_), size_u(size()));
        }


        auto uu()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "uu", get(vertex_index, graph_), size_u(size()));
        }


        auto uu() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "uu", get(vertex_index, graph_), size_u(size()));
        }


        auto C()
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "C", get(vertex_index, graph_), size_C(size()));
        }


        auto C() const
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "C", get(vertex_index, graph_), size_C(size()));
        }


        auto D()
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "D", get(vertex_index, graph_), size_D(size()));
        }


        auto D() const
        {
            return detail::makeJsonQpPropertyMap<OcpVertexDescriptor, DynamicMatrix<Kernel>>(
                json_.at("nodes"), "D", get(vertex_index, graph_), size_D(size()));
        }


        auto ld()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ld", get(vertex_index, graph_), size_d(size()));
        }


        auto ld() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ld", get(vertex_index, graph_), size_d(size()));
        }


        auto ud()
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ud", get(vertex_index, graph_), size_d(size()));
        }


        auto ud() const
        {
            return detail::makeJsonQpPropertyMap<vertex_descriptor, DynamicVector<Kernel>>(
                json_.at("nodes"), "ud", get(vertex_index, graph_), size_d(size()));
        }


        auto A()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "A", get(edge_index, graph_), size_A(size(), graph_));
        }


        auto A() const
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "A", get(edge_index, graph_), size_A(size(), graph_));
        }


        auto B()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "B", get(edge_index, graph_), size_B(size(), graph_));
        }


        auto B() const
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicMatrix<Kernel>>(
                json_["edges"], "B", get(edge_index, graph_), size_B(size(), graph_));
        }


        auto b()
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicVector<Kernel>>(
                json_["edges"], "b", get(edge_index, graph_), size_b(size(), graph_));
        }


        auto b() const
        {
            return detail::makeJsonQpPropertyMap<edge_descriptor, DynamicVector<Kernel>>(
                json_["edges"], "b", get(edge_index, graph_), size_b(size(), graph_));
        }


    private:
        // OCP graph
        OcpGraph graph_;

        // OCP node sizes
        std::vector<OcpSize> size_;

        // OCP JSON
        tmpc::json json_;
    };


    template <typename Kernel>
    struct KernelOf<JsonQp<Kernel>>
    {
        using type = Kernel;
    };


    template <typename Kernel>
    struct RealOf<JsonQp<Kernel>>
    {
        using type = typename Kernel::Real;
    };
}
