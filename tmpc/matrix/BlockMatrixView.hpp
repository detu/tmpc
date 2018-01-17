#pragma once

#include <boost/range/iterator_range.hpp>

#include <vector>
#include <stdexcept>


namespace tmpc
{
    template <typename Kernel, typename MT>
    class BlockMatrixView
    {
    public:
        using size_t = typename Kernel::size_t;


        BlockMatrixView(MT& m, std::initializer_list<size_t> size_m, std::initializer_list<size_t> size_n)
        :   BlockMatrixView(m, boost::make_iterator_range(size_m), boost::make_iterator_range(size_n))
        {            
        }

        
        template <typename IteratorRange1, typename IteratorRange2>
        BlockMatrixView(MT& m, IteratorRange1 const& size_m, IteratorRange2 const& size_n)
        :   m_ {m}
        ,   rowIndex_(size_m.size() + 1)
        ,   colIndex_(size_n.size() + 1)
        {
            auto si = size_m.begin();
            auto sj = size_n.begin();

            rowIndex_[0] = 0;
            for (std::size_t i = 1; i < rowIndex_.size(); ++i)
                rowIndex_[i] = rowIndex_[i - 1] + *si++;

            colIndex_[0] = 0;
            for (std::size_t i = 1; i < colIndex_.size(); ++i)
                colIndex_[i] = colIndex_[i - 1] + *sj++;

            if (rowIndex_.back() != rows(m))
                throw std::invalid_argument("Sum of block row count must equal to matrix row count");

            if (colIndex_.back() != columns(m))
                throw std::invalid_argument("Sum of block column count must equal to matrix column count");
        }

        
        friend auto rows(BlockMatrixView const& block_view)
        {
            return block_view.rowIndex_.size() - 1;
        }


        friend auto columns(BlockMatrixView const& block_view)
        {
            return block_view.colIndex_.size() - 1;
        }


        decltype(auto) operator()(size_t i, size_t j) const
        {
            assert(i + 1 < rowIndex_.size() && j + 1 < colIndex_.size());
            return submatrix(m_, rowIndex_[i], colIndex_[j], rowIndex_[i + 1] - rowIndex_[i], colIndex_[j + 1] - colIndex_[j]);
        }


    private:
        MT& m_;
        std::vector<size_t> rowIndex_;
        std::vector<size_t> colIndex_;
    };
}