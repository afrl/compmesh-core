#pragma once

#include "composite_mesher/FeatureIdentification.h"
#include "composite_mesher/GlobalOptions.h"
#include "composite_mesher/Tetrahedralization.h"
#include "mesh/algorithms/IO.h"
#include "simple_math/Utility.h"

namespace composite_mesher
{
struct UpdatedEdgeInfo
{
    enum Type
    {
        NEW_EDGE,
        SUBDIVIDED_EDGE
    };

    Type update_type;
    mesh::LocalIndex orig_edge_id;
    mesh::LocalIndex new_edge_id;
    mesh::LocalIndex new_node_id;
};

struct UpdatedTriangleInfo
{
    mesh::LocalIndex orig_tri_id;
    std::vector<mesh::LocalIndex> new_tri_ids;
};

struct UpdatedTetrahedronInfo
{
    mesh::LocalIndex orig_tet_id;
    std::vector<mesh::LocalIndex> new_tet_ids;
};

/// @brief Stores information about how a tetrahedralization is updated during remeshing
struct UpdateResult
{
    std::vector<UpdatedEdgeInfo> updated_edges;
    std::vector<UpdatedTriangleInfo> updated_triangles;
    std::vector<UpdatedTetrahedronInfo> updated_tetrahedra;
};

/// @brief Splits all of the given edges in the tetrahedralization, subdividing edges, triangle
/// faces, and tetrahedra as needed.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         Tetrahedralization& t);

/// @brief Splits all of the given edges in the tetrahedralization, subdividing edges, triangle
/// faces, and tetrahedra as needed but uses the given existing vertices for the split edge.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         const std::vector<mesh::LocalIndex>& existing_vertices_for_split_edges,
                         Tetrahedralization& t);
}