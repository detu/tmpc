#pragma once

#include <vector>
#include <iterator>

namespace tmpc {
namespace detail {

template <typename Base_>
class MultiStage : public Base_
{
public:
    typedef Base_ Base;

    template <typename InputIterator>
    MultiStage(InputIterator so_first, InputIterator so_last)
    :   Base(so_first, so_last)
    {
        stage_.reserve(std::distance(so_first, so_last));
        for (InputIterator current = so_first; current != so_last; ++current)
            stage_.push_back(Base::createStage(current, so_first, so_lase));
    }

    Stage& operator[](std::size_t i)
    {
        return stage_.at(i);
    }

    Stage const& operator[](std::size_t i) const
    {
        return stage_.at(i);
    }

    std::size_t size() const
    {
        return stage_.size();
    }

    typedef typename std::vector<Stage>::iterator iterator;
    typedef typename std::vector<Stage>::const_iterator const_iterator;
    typedef typename std::vector<Stage>::reference reference;
    typedef typename std::vector<Stage>::const_reference const_reference;

    iterator begin()
    {
        return stage_.begin();
    }

    iterator end()
    {
        return stage_.end();
    }

    const_iterator begin() const
    {
        return stage_.begin();
    }

    const_iterator end() const
    {
        return stage_.end();
    }

    reference front()
    {
        return stage_.front();
    }

    reference back()
    {
        return stage_.back();
    }

    const_reference front() const
    {
        return stage_.front();
    }

    const_reference back() const
    {
        return stage_.back();
    }

private:
    std::vector<Stage> stage_;
};

}   // namespace detail
}   // namespace tmpc