/// \brief Some useful math functions.

#pragma once

#include <blaze/Math.h>

#include <boost/throw_exception.hpp>

#include <limits>
#include <cmath>
#include <stdexcept>


namespace tmpc
{
    using std::cos;
    using std::sin;
    using std::pow;
    using std::exp;
    using std::abs;


    template <typename T>
    inline T constexpr inf()
    {
        return std::numeric_limits<T>::infinity();
    }


    template <typename T>
    inline T constexpr sNaN()
    {
        return std::numeric_limits<T>::signaling_NaN();
    }


    namespace detail
    {
        template <typename T>
        inline size_t maxRowLength(std::initializer_list<std::initializer_list<T>> const& l)
        {
            size_t m = 0;

            for (auto const& row : l)
                m = std::max(m, size(row));

            return m;
        }


        template <typename T>
        class NoResize;


        template <typename ET, bool SO>
        class NoResize<blaze::DynamicMatrix<ET, SO>>
        {
        public:
            using MatrixType = blaze::DynamicMatrix<ET, SO>;


            explicit NoResize(MatrixType& m)
            :   m_(m)
            {
            }


            template <typename MT1, bool SO1>
            MatrixType& operator=(blaze::Matrix<MT1, SO1> const& rhs)
            {
                if (rows(rhs) != rows(m_) || columns(rhs) != columns(m_))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a matrix assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return m_ = rhs;
            }


            template <typename T>
		    MatrixType& operator=(std::initializer_list<std::initializer_list<T>> rhs)
            {
                if (!(size(rhs) == rows(m_) && (size(rhs) == 0 || maxRowLength(rhs) == columns(m_))))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a matrix assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return m_ = rhs;
            }


            template <typename T>
            MatrixType& operator=(T const& rhs)
            {
                return m_ = rhs;
            }


        private:
            MatrixType& m_;
        };


        template <typename ET, bool SO>
        class NoResize<blaze::SymmetricMatrix<blaze::DynamicMatrix<ET, SO>>>
        {
        public:
            using MatrixType = blaze::SymmetricMatrix<blaze::DynamicMatrix<ET, SO>>;


            explicit NoResize(MatrixType& m)
            :   m_(m)
            {
            }


            template <typename MT1, bool SO1>
            MatrixType& operator=(blaze::Matrix<MT1, SO1> const& rhs)
            {
                if (rows(rhs) != rows(m_) || columns(rhs) != columns(m_))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a matrix assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return m_ = rhs;
            }


            template <typename T>
		    MatrixType& operator=(std::initializer_list<std::initializer_list<T>> rhs)
            {
                if (!(size(rhs) == rows(m_) && (size(rhs) == 0 || maxRowLength(rhs) == columns(m_))))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a matrix assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return m_ = rhs;
            }


            template <typename T>
            MatrixType& operator=(T const& rhs)
            {
                return m_ = rhs;
            }


        private:
            MatrixType& m_;
        };


        template <typename ET, bool TF>
        class NoResize<blaze::DynamicVector<ET, TF>>
        {
        public:
            using VectorType = blaze::DynamicVector<ET, TF>;


            explicit NoResize(VectorType& v)
            :   v_(v)
            {
            }


            template <typename VT1, bool TF1>
            VectorType& operator=(blaze::Vector<VT1, TF1> const& rhs)
            {
                if (size(rhs) != size(v_))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a vector assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return v_ = rhs;
            }


            template <typename T>
		    VectorType& operator=(std::initializer_list<T> rhs)
            {
                if (size(rhs) != size(v_))
                    BOOST_THROW_EXCEPTION(std::invalid_argument("Right-hind side of a vector assignment has different size, "
                        "but resizing of the left-hand size is not allowed"));

                return v_ = rhs;
            }


            template <typename T>
            VectorType& operator=(T const& rhs)
            {
                return v_ = rhs;
            }


        private:
            VectorType& v_;
        };
    }


    // TODO: uncomment the following lines when 
    // we finally get rid of Eigen wrappers with the same names.
    //
    // using blaze::StaticMatrix;
    // using blaze::StaticVector;
    // using blaze::DynamicMatrix;
    // using blaze::DynamicVector;
    // using blaze::columnMajor;
    // using blaze::rowMajor;
    // using blaze::columnVector;
    // using blaze::rowVector;
}


namespace blaze
{
    /// @brief No-resize adaptor for vector assignment
    ///
    /// See https://bitbucket.org/blaze-lib/blaze/issues/93/matrix-vector-assignment-without-resizing
    template <typename ET, bool TF>
    inline auto noresize(DynamicVector<ET, TF>& v)
    {
        return tmpc::detail::NoResize<DynamicVector<ET, TF>>(v);
    }


    /// @brief No-resize adaptor for matrix assignment
    ///
    /// See https://bitbucket.org/blaze-lib/blaze/issues/93/matrix-vector-assignment-without-resizing
    template <typename ET, bool SO>
    inline auto noresize(DynamicMatrix<ET, SO>& m)
    {
        return tmpc::detail::NoResize<DynamicMatrix<ET, SO>>(m);
    }


    /// @brief No-resize adaptor for symmetric matrix assignment
    ///
    /// See https://bitbucket.org/blaze-lib/blaze/issues/93/matrix-vector-assignment-without-resizing
    template <typename ET, bool SO>
    inline auto noresize(SymmetricMatrix<DynamicMatrix<ET, SO>>& m)
    {
        return tmpc::detail::NoResize<SymmetricMatrix<DynamicMatrix<ET, SO>>>(m);
    }


    template <typename VT, bool TF>
    inline decltype(auto) scalar(Vector<VT, TF> const& v)
    {
        if (size(v) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a vector of size 1 expected"));

        return (~v)[0];
    }


    template <typename VT, bool TF>
    inline decltype(auto) scalar(Vector<VT, TF>& v)
    {
        if (size(v) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a vector of size 1 expected"));

        return (~v)[0];
    }


    template <typename VT, bool TF>
    inline decltype(auto) scalar(Vector<VT, TF>&& v)
    {
        if (size(v) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a vector of size 1 expected"));

        return (~v)[0];
    }


    template <typename MT, bool SO>
    inline decltype(auto) scalar(Matrix<MT, SO> const& m)
    {
        if (rows(m) != 1 || columns(m) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a 1x1 matrix expected"));

        return (~m)(0, 0);
    }


    template <typename MT, bool SO>
    inline decltype(auto) scalar(Matrix<MT, SO>& m)
    {
        if (rows(m) != 1 || columns(m) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a 1x1 matrix expected"));

        return (~m)(0, 0);
    }


    template <typename MT, bool SO>
    inline decltype(auto) scalar(Matrix<MT, SO>&& m)
    {
        if (rows(m) != 1 || columns(m) != 1)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Invalid size of the scalar() argument: a 1x1 matrix expected"));

        return (~m)(0, 0);
    }
}