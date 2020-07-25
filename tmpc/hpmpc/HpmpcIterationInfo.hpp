#pragma once

#include <array>


namespace tmpc
{
    /// @brief HPMPC iteration info.
    ///
    /// HPMPC returns 5 Real numbers per iteration:
    /// step length for predictor and corrector, 
    /// centering parameter, duality measure for predictor and corrector.
    template <typename Real>
    struct HpmpcIterationInfo
    {
        /// @brief Step length for predictor
        Real stepLengthPredictor;

        /// @brief Step length for corrector
        Real stepLengthCorrector;

        /// @brief Centering parameter
        Real centeringParameter;

        /// @brief Duality measure for predictor
        Real dualityMeasurePredictor;

        /// @brief Duality measure for corrector
        Real dualityMeasureCorrector;
    };
}