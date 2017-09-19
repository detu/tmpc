#pragma once

namespace tmpc
{
    enum TransposeFlag
    {
        columnVector,
        rowVector
    };

    TransposeFlag constexpr defaultTransposeFlag = columnVector;
}