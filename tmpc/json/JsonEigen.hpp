#pragma once

#include "Json.hpp"

#include <Eigen/Dense>

#include <algorithm>


namespace Eigen
{
    template <typename T, int Rows, int Cols, int Options>
    void to_json(tmpc::json& jsn, Matrix<T, Rows, Cols, Options> const& m) 
    {
        jsn = tmpc::json::array();

        /*
        for (size_t i = 0; i < m.rows(); ++i)
        {
            tmpc::json row = tmpc::json::array();

            for (size_t j = 0; j < columns(m); ++j)
                row.push_back((~m)(i, j));

            jsn.push_back(row);
        }
        */
    }


    template <typename VT>
    void to_json(tmpc::json& jsn, MatrixBase<VT> const& v) 
    {
        jsn = tmpc::json::array();
        //std::copy(v.begin(), v.end(), std::back_inserter(jsn));
    }


    template <typename VT>
    void to_json(tmpc::json& jsn, PlainObjectBase<VT> const& v) 
    {
        jsn = tmpc::json::array();
        //std::copy(v.begin(), v.end(), std::back_inserter(jsn));
    }


    // TODO: change from_json() function for matrices such that they expect
    // json array of arrays as input.
    template <typename VT>
    void from_json(tmpc::json const& j, PlainObjectBase<VT>& v)
    {
        using Scalar = typename VT::Scalar;

        v.resize(j.size());
        std::transform(j.begin(), j.end(), v.data(), [] (nlohmann::json const& val) -> Scalar
        {
            if (val.is_string())
            {
                static constexpr Scalar inf = std::numeric_limits<Scalar>::infinity();

                if (val == "-inf")
                    return -inf;
                else if (val == "inf")
                    return inf;
                else
                {
                    std::ostringstream msg;
                    msg << "Invalid value in config file: a double number, a \"inf\" or a \"-inf\" was expected, but \"" << val << "\" found.";
                    throw std::invalid_argument(msg.str());
                }
            }
            else
                return val;
        });
    }

    template <typename T, int Rows, int Cols, int Options>
    void from_json(tmpc::json const& j, Matrix<T, Rows, Cols, Options>& v)
    {
        using Scalar = T;
        
        if (v.size() != j.size())
            throw std::invalid_argument("Number of elements in json value does not match matrix size");

        int ii = 0, jj = 0;
        for (auto const& val : j)
        {
            T numeric_val = std::numeric_limits<Scalar>::signaling_NaN();

            if (val.is_string())
            {
                static constexpr Scalar inf = std::numeric_limits<Scalar>::infinity();

                if (val == "-inf")
                    numeric_val = -inf;
                else if (val == "inf")
                    numeric_val = inf;
                else
                {
                    std::ostringstream msg;
                    msg << "Invalid value in config file: a number, a \"inf\" or a \"-inf\" was expected, but \"" << val << "\" found.";
                    throw std::invalid_argument(msg.str());
                }
            }
            else
                numeric_val = val;

            v(ii, jj) = numeric_val;

            if (++jj >= v.cols())
            {
                jj = 0;                
                ++ii;
            }
        }
    }
}