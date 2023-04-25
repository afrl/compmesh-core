#pragma once

#include "Seams.h"
#include "Tetrahedralization.h"
#include "Triangulation.h"

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

// Split edges for Tetrahedralizations

/// @brief Splits all of the given edges in the tetrahedralization, subdividing edges, triangle
/// faces, and tetrahedra as needed.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         Tetrahedralization& t);

/// @brief Splits all of the given edges in the tetrahedralization, subdividing edges, triangle
/// faces, and tetrahedra as needed but uses the given existing vertices for the split edge.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         const std::vector<mesh::LocalIndex>& existing_vertices_for_split_edges,
                         Tetrahedralization& t);

// Split edges for Triangulations

/// @brief Splits all of the given edges in the triangulation, subdividing edges and triangles
/// as needed.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split, Triangulation& t);

/// @brief Splits all of the given edges in the triangulation, subdividing edges and triangles
/// as needed but uses the given existing vertices for the split edge.
UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         const std::vector<mesh::LocalIndex>& existing_vertices_for_split_edges,
                         Triangulation& t);

/// @brief Split all edges for which the two vertices lie on any seam but the edge does not lie on a
///        seam.
///
/// @note Is called by remove_interpenetrations
///
/// When a Triangulation is optimize, tangential smoothing creates problems if constrained edges
/// exist.  So this is important to call before optimization.
///
/// For example, the edges e1 and e2 will be split:
///             seam edge
///             o-------o
///            / \     /
/// seam edge /   \e1 /e2
///          /     \ /
///         o-------o
///         seam edge
///
UpdateResult split_constrained_edges(Triangulation& t, SeamInfo& seams);
}