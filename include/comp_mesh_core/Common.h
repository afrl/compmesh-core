#pragma once

#include "mesh/Defines.h"

namespace composite_mesher
{
using ConnIndex = uint8_t;
using Index     = uint32_t;
using Marker    = uint32_t;
using FacetID   = Marker;

/// @brief A fast, safe hash for linear 2-node edges
inline const uint64_t hash_linear_edge(const Index& v1, const Index& v2)
{
    if (v1 < v2)
    {
        uint64_t ret = 0;
        ret |= v1;
        ret |= static_cast<int64_t>(v2) << 32;
        return ret;
    }
    else
    {
        uint64_t ret = 0;
        ret |= v2;
        ret |= static_cast<int64_t>(v1) << 32;
        return ret;
    }
};
}