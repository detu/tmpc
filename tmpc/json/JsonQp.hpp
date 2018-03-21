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
}
