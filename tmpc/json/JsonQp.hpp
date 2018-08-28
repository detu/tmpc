/*
@brief Support convertion between OCP QP and JSON.
*/

#pragma once

#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/Matrix.hpp>

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
    }


    template <typename Kernel>
    class JsonQp
    {
    public:
        using vertex_descriptor = OcpVertexDescriptor;
        using edge_descriptor = OcpEdgeDescriptor;


        JsonQp()
        {
        }


        auto& graph()
        {
            return graph_;
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto const& json() const
        {
            return json_;
        }


        template <typename KeyType, typename ValueType>
        struct PropertyMap
        {
            using key_type = KeyType; 
            using value_type = ValueType;
            using category = boost::read_write_property_map_tag;

            PropertyMap(JsonQp& qp, std::string const& name) 
            :   qp_{qp}
            ,   name_{name}
            {
            }


            JsonQp& qp_;
            std::string const name_;
        };


        auto Q()
        {
            return PropertyMap<vertex_descriptor, DynamicMatrix<Kernel>> {*this, "Q"};
        }


        auto R()
        {
            return PropertyMap<vertex_descriptor, DynamicMatrix<Kernel>> {*this, "R"};
        }


        auto S()
        {
            return PropertyMap<vertex_descriptor, DynamicMatrix<Kernel>> {*this, "S"};
        }


        auto q()
        {
            return PropertyMap<vertex_descriptor, DynamicVector<Kernel>> {*this, "q"};
        }


        auto r()
        {
            return PropertyMap<vertex_descriptor, DynamicVector<Kernel>> {*this, "r"};
        }


        auto A()
        {
            return PropertyMap<edge_descriptor, DynamicMatrix<Kernel>> {*this, "A"};
        }


        auto B()
        {
            return PropertyMap<edge_descriptor, DynamicMatrix<Kernel>> {*this, "B"};
        }


        auto b()
        {
            return PropertyMap<edge_descriptor, DynamicVector<Kernel>> {*this, "b"};
        }


        template <typename KeyType, typename ValueType>
        friend void put(PropertyMap<KeyType, ValueType> const& m, KeyType const& key, ValueType const& val)
        {
            using traits = detail::DescriptorTraits<KeyType>;
            auto const top_key = traits::dictionaryKey();
            auto const index = get(typename traits::index_property_t(), m.qp_.graph_);
            m.qp_.json_[top_key][get(index, key)][m.name_] = val;
        }


        template <typename KeyType, typename ValueType>
        friend ValueType get(PropertyMap<KeyType, ValueType> const& m, KeyType const& key)
        {
            using traits = detail::DescriptorTraits<KeyType>;
            auto const top_key = traits::dictionaryKey();
            auto const index = get(typename traits::index_property_t(), m.qp_.graph_);
            
            return m.qp_.json_[top_key][get(index, key)][m.name_];
        }


    private:
        OcpGraph graph_;
        ::tmpc::json json_;
    };
}
