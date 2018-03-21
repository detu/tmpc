#pragma once

#include "Json.hpp"

#include <blaze/Blaze.h>


namespace blaze
{
    template <typename MT, bool SO>
    void to_json(tmpc::json& jsn, Matrix<MT, SO> const& m) 
    {
        jsn = tmpc::json::array();

        for (size_t i = 0; i < rows(m); ++i)
            for (size_t j = 0; j < columns(m); ++j)
                jsn.push_back((~m)(i, j));
    }


    template <typename VT, bool SO>
    void to_json(tmpc::json& jsn, Vector<VT, SO> const& v) 
    {
        jsn = tmpc::json::array();
        std::copy((~v).begin(), (~v).end(), std::back_inserter(jsn));
    }


    template <typename MT, bool SO>
    void from_json(tmpc::json const& j, Matrix<MT, SO>& v)
    {
        using Scalar = typename MT::ElementType;
        
        if (rows(v) * columns(v) != j.size())
            throw std::invalid_argument("Number of elements in json value does not match matrix size");

        int ii = 0, jj = 0;
        for (auto const& val : j)
        {
            auto numeric_val = std::numeric_limits<Scalar>::signaling_NaN();

            if (val.is_string())
            {
                static auto constexpr inf = std::numeric_limits<Scalar>::infinity();

                if (val == "-inf")
                    numeric_val = -inf;
                else if (val == "inf")
                    numeric_val = inf;
                else
                {
                    std::ostringstream msg;
                    msg << "Invalid floating point value in json file: a number, a \"inf\" or a \"-inf\" was expected, but \"" << val << "\" found.";
                    throw std::invalid_argument(msg.str());
                }
            }
            else
                numeric_val = val;

            (~v)(ii, jj) = numeric_val;

            if (++jj >= columns(v))
            {
                jj = 0;                
                ++ii;
            }
        }
    }
    

    template <typename VT, bool SO>
    void from_json(tmpc::json const& j, Vector<VT, SO>& v) 
    {
        using Scalar = typename VT::ElementType;

        v.resize(j.size());
        std::transform(j.begin(), j.end(), (~v).begin(), [] (tmpc::json const& val) -> Scalar
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
                    msg << "Invalid floating point value in json file: a number, a \"inf\" or a \"-inf\" was expected, but \"" << val << "\" found.";
                    throw std::invalid_argument(msg.str());
                }
            }
            else
                return val;
        });
    }
}