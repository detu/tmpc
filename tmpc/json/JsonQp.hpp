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


    /*
    template <typename QP>
    void from_json(json const& j, OcpQpBase<QP>& qp) 
    {
        // Not implemented        
    }
    */

    template <typename Kernel>
    void from_json(json const& j, OcpQp<Kernel> qp)
    {
        OcpSize sz {j["q"].size(), j["r"].size(), j["ld"].size(), j["zl"].size()};
        qp = OcpQp<Kernel> {sz};

        DynamicMatrix<Kernel> Q {sz.nx(), sz.nx()};
        from_json(j["Q"], Q);
        qp.Q(Q);
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
