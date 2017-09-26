#pragma once

namespace tmpc
{
    enum StorageOrder
    {
        rowMajor,
        columnMajor
    };

    static StorageOrder constexpr defaultStorageOrder = columnMajor;
}