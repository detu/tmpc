#pragma once

#include <tmpc/ocp/OcpPoint.hpp>

#include <boost/range/iterator_range_core.hpp>

#include <vector>


namespace tmpc
{
    template <typename QpWorkspace>
    class SqpWorkspace
    {
    public:
        template <typename IteratorRange>
        SqpWorkspace(IteratorRange const& ocp_size)
        :   qpWorkspace_ {ocp_size}
        ,   point_ {ocp_size.size()}
        {
        }

        decltype(auto) qp()
        {
            return qpWorkspace_.qp();
        }

        decltype(auto) qp() const
        {
            return qpWorkspace_.qp();
        }

        decltype(auto) step()
        {
            qpWorkspace_.solve();
            return qpWorkspace_.solution();
        }

    private:
        QpWorkspace qpWorkspace_;
        std::vector<OcpPoint> point_;
    };
}