#pragma once

#include "mesh/Defines.h"

namespace composite_mesher
{
using ConnIndex = mesh::ConnIndex;
using Index     = mesh::LocalIndex;
using Marker    = uint32_t;
using FacetID   = Marker;
using SeamID   = Marker;

constexpr Index NOT_SET_Index   = std::numeric_limits<Index>::max();
constexpr Marker NOT_SET_Marker = std::numeric_limits<Marker>::max();
constexpr FacetID NOT_SET_FacetID = std::numeric_limits<FacetID>::max();
constexpr SeamID NOT_SET_SeamID   = std::numeric_limits<SeamID>::max();

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