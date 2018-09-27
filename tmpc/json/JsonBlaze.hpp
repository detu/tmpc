#pragma once

#include "Json.hpp"

#include <blaze/Blaze.h>


namespace blaze
{
    template <typename MT, bool SO>
    void to_json(tmpc::json& jsn, Matrix<MT, SO> const& m) 
    {
        jsn = tmpc::json::array();

        // Encoding empty matrices just as empty json arrays `[]`.
        // Although the dimension information can be preserved for Nx0 matrices,
        // for 0xN matrices it can not.
        if (!isEmpty(m))
        {
            for (size_t i = 0; i < rows(m); ++i)
            {
                tmpc::json row = tmpc::json::array();

                for (size_t j = 0; j < columns(m); ++j)
                    row.push_back((~m)(i, j));

                jsn.push_back(row);
            }
        }
    }


    template <typename VT, bool SO>
    void to_json(tmpc::json& jsn, Vector<VT, SO> const& v) 
    {
        jsn = tmpc::json::array();
        std::copy((~v).begin(), (~v).end(), std::back_inserter(jsn));
    }


    template <typename MT, bool SO>
    void from_json(tmpc::json const& j, Matrix<MT, SO>& m)
    {
        using Scalar = typename MT::ElementType;

        if (j.is_number())
        {
            // Interpret number as a 1x1 matrix
            (~m).resize(1, 1);
            (~m)(0, 0) = j;
        }
        else
        {
            size_t const M = j.size();
            size_t const N = M > 0 ? j[0].size() : 0;

            (~m).resize(M, N);
            
            for (size_t ii = 0; ii < M; ++ii)
            {
                if (j[ii].size() != N)
                    throw std::invalid_argument("Non-equal number of columns in json matrix");

                size_t jj = 0;
                for (auto val : j[ii])
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

                    (~m)(ii, jj) = numeric_val;
                    ++jj;
                }
            }
        }
    }
    

    template <typename VT, bool SO>
    void from_json(tmpc::json const& j, Vector<VT, SO>& v) 
    {
        using Scalar = typename VT::ElementType;

        /*
        if (j.is_number())
        {
            // Interpret numbers as 1-dimensional vectors
            (~v).resize(1);
            (~v)[0] = j;
        }
        */

        (~v).resize(j.size());
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