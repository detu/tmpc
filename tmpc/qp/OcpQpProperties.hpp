#pragma once

#include <tmpc/ocp/OcpGraph.hpp>

#include <blaze/Math.h>

namespace tmpc
{
    /// @brief Set A, B and b to represent a linearized shooting equality
    /// of the form \Delta x_{k+1} = \frac{\dif f}{\dif x}(x_k,u_k)\, \Delta x_{k} 
    ///	+ \frac{\dif f}{\dif u}(x_k,u_k)\, \Delta u_{k} + f(x_{k},u_{k}) - x_{k+1}
    ///
    /// @param qp QP object
    /// @param e QP edge for which to set the properties
    /// @param f system state at target node according to system dynamics
    /// @param Jx derivative df/dx
    /// @param Ju derivative df/du
    /// @param map_x mapping from a vertex to a linearization point in state space
    template <typename QP, typename Vector1, typename Matrix1, typename Matrix2, typename MapX>
    inline void linearizedShootingEquality(QP& qp, OcpEdgeDescriptor e, 
        Vector1 const& f, Matrix1 const& Jx, Matrix2 const& Ju, MapX map_x)
    {
        put(qp.A(), e, Jx);
        put(qp.B(), e, Ju);
        put(qp.b(), e, f - get(map_x, target(e, qp.graph())));
    }


    /// @brief Set upper and lower input bounds relative to a point.
    ///
    /// TODO: think about making u(), x() property maps of a new SqpWorkspace class
    /// (see https://gitlab.syscop.de/mikhail.katliar/tmpc/issues/45)
    template <typename QP, typename MapU, typename Vector2, typename Vector3>
    inline void relativeInputBounds(QP& qp, OcpVertexDescriptor v, MapU u, Vector2 const& lu, Vector3 const& uu)
    {
        put(qp.lu(), v, lu - get(u, v));
        put(qp.uu(), v, uu - get(u, v));
    }


    /// @brief Set upper and lower state bounds relative to a point.
    template <typename QP, typename MapX, typename Vector2, typename Vector3>
    inline void relativeStateBounds(QP& qp, OcpVertexDescriptor v, MapX x, Vector2 const& lx, Vector3 const& ux)
    {
        put(qp.lx(), v, lx - get(x, v));
        put(qp.ux(), v, ux - get(x, v));
    }


    /// @brief Set the Gauss-Newton approximation of the hessian Hessian and the gradient.
    template <typename QP, typename ResidualVector, typename CMatrix, typename DMatrix>
    inline void gaussNewtonCostApproximation(QP& qp, OcpVertexDescriptor v, ResidualVector const& res, CMatrix const& C, DMatrix const& D)
    {
        // H = G^T G
        //   = [Q S^T
        //      S R]
        //

        put(qp.Q(), v, trans(C) * C);
        put(qp.R(), v, trans(D) * D);
        put(qp.S(), v, trans(D) * C);

        // g = 2 * (y_bar - y_hat)^T * W * G
        // g = [q; r]
        put(qp.q(), v, trans(C) * res);
        put(qp.r(), v, trans(D) * res);
    }


    ///@brief Add Levenberg-Marquardt term to Q and R.
    template <typename QP, typename Real>
    inline void addLevenbergMarquardt(QP& qp, OcpVertexDescriptor v, Real levenberg_marquardt)
    {
        put(qp.Q(), v, get(qp.Q(), v) + levenberg_marquardt * blaze::IdentityMatrix<Real>(get(qp.size(), v).nx()));
        put(qp.R(), v, get(qp.R(), v) + levenberg_marquardt * blaze::IdentityMatrix<Real>(get(qp.size(), v).nu()));
    }
}