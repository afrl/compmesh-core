#pragma once

#include "Common.h"
#include "Triangulation.h"

namespace composite_mesher
{

// Forward declare (see below)
struct StaticVertexMapAccessor;

/// @brief Used to store a vertex-to-edge, vertex-to-triangle, and vertex-to-tetrahedra maps
/// Note: this struct is designed to lookup performance, but it is expensive to update as a result
/// (hence the struct name).  If you plan to remove or renumber vertices, then you will want to make
/// a different structure for those cases.
struct StaticVertexMap
{
    std::vector<Index> id_list;
    std::vector<std::size_t> offsets;

    StaticVertexMapAccessor get_indices(const Index& vertex_id) const;
};

/// @brief Helper for StaticVertexMap to enable range-based loops with StaticVertexMap::get_indices
struct StaticVertexMapAccessor
{
    auto begin() const
    {
        return map.id_list.begin() + map.offsets[vertex_id];
    }
    auto end() const
    {
        return map.id_list.begin() + map.offsets[vertex_id + 1];
    }

    const StaticVertexMap& map;
    const Index vertex_id;
};

StaticVertexMap build_vert_to_tri_map(const Triangulation& t);
}