#pragma once


namespace tmpc
{
    /// @brief Defines the OCP KKT value concept.
    ///
    template <typename KktValue>
    concept OcpKktValue = requires(KktValue kkt, OcpVertex v, OcpEdge e)
    {
        typename KktValue::Real;

        kkt.graph();
        
        kkt.gx(v);
        kkt.gu(v);
        kkt.c(e);
    };
}