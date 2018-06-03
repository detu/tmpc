/*
@brief Support convertion between OCP QP and JSON.
*/

#pragma once

#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/OcpQp.hpp>

#include "Json.hpp"


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
    void from_json(json const& j, OcpQpBase<QP>& qp)
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
}
