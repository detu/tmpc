#pragma once

#include <iterator>

namespace tmpc {
namespace detail {

template <typename T, typename Container = std::vector<T>>
class NonResizableCollection
{
public:
    typedef T Element;

    Element& operator[](std::size_t i)
    {
        return elements_.at(i);
    }

    Element const& operator[](std::size_t i) const
    {
        return elements_.at(i);
    }

    std::size_t size() const
    {
        return elements_.size();
    }

    typedef typename Container::iterator iterator;
    typedef typename Container::const_iterator const_iterator;
    typedef typename Container::reference reference;
    typedef typename Container::const_reference const_reference;

    iterator begin()
    {
        return elements_.begin();
    }

    iterator end()
    {
        return elements_.end();
    }

    const_iterator begin() const
    {
        return elements_.begin();
    }

    const_iterator end() const
    {
        return elements_.end();
    }

    reference front()
    {
        return elements_.front();
    }

    reference back()
    {
        return elements_.back();
    }

    const_reference front() const
    {
        return elements_.front();
    }

    const_reference back() const
    {
        return elements_.back();
    }

protected:
    Container& elements()
    {
        return elements_;
    }

    Container const& elements() const
    {
        return elements_;
    }

private:
    Container elements_;
};

template <typename T, typename Container = std::vector<T>>
class NonResizableClassAdapter
{
public:
    typedef T Element;

    Element& operator[](std::size_t i)
    {
        return elements_.at(i);
    }

    Element const& operator[](std::size_t i) const
    {
        return elements_.at(i);
    }

    std::size_t size() const
    {
        return elements_.size();
    }

    typedef typename Container::iterator iterator;
    typedef typename Container::const_iterator const_iterator;
    typedef typename Container::reference reference;
    typedef typename Container::const_reference const_reference;

    iterator begin()
    {
        return elements_.begin();
    }

    iterator end()
    {
        return elements_.end();
    }

    const_iterator begin() const
    {
        return elements_.begin();
    }

    const_iterator end() const
    {
        return elements_.end();
    }

    reference front()
    {
        return elements_.front();
    }

    reference back()
    {
        return elements_.back();
    }

    const_reference front() const
    {
        return elements_.front();
    }

    const_reference back() const
    {
        return elements_.back();
    }

protected:
    Container& elements()
    {
        return elements_;
    }

    Container const& elements() const
    {
        return elements_;
    }

private:
    Container elements_;
};

}   // namespace detail
}   // namespace tmpc