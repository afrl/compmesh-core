#pragma once

#include "comp_mesh_core/Triangulation.h"

namespace composite_mesher
{
struct CollapseResult
{
    std::vector<Index> removed_vertices;
    std::vector<Index> removed_edges;
    std::vector<Index> removed_triangles;
};

/// @brief Collapses the given edges in the triangulation, merging the two vertices into one (which
/// will be at the midpoint).
CollapseResult collapse_edges(const std::vector<Index>& edges_to_collapse, Triangulation& t);

}