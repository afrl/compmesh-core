#pragma once
#pragma once

#include "Eigen/Eigen"
#include "Tetrahedralization.h"

#include <array>
#include <set>
#include <unordered_map>

namespace composite_mesher
{

/// @brief A new extended marching tetrahedra algorithm to support isosurfaces that exactly
/// intersect vertices along the boundary of the tetrahedra grid.
///
/// Marching tetrahedra is an established algorithm to find an isosurface given a continuous
/// function evaluated on a tetrahedra grid.
/// @param tet_vertices - coordinates of vertices organized as Nx3 matrix, where N is the number of
///        vertices
/// @param tet_connectivity - connectivity of tetrahedra as Tx4 matrix (where T is the number of
///        tetrahedra) following numbering:
///         Front face        Top down
///             o 3              o 0
///            /|\              / \ 
///           / | \            / o \ 
///          o--o--o          o-----o
///         1   0   2         1     2
/// @param isovalues - continuous function value at each given vertex (length must match
///        tet_vertices.rows())
/// @out tri_vertices - vertices for the output isosurface (may contain points that appear in
///      tet_vertices)
/// @out tri_connectivity  - connectivity of triangles along output isosurface (numbering is
///      counter-clock and normal points along positive direction of continuous function field)
/// @out tri_vertex_to_tet_vertex_map - a map from the output vertex index to the same vertex in
///      tet_vertices (if there is one, otherwise, the value will be -1 for new vertices)
/// @out tri_vertex_to_tet_split_edge_map - a map from the vertex to the split edge on which it
///      was created (if there is one, otherwise, the row will be -1 for existing vertices)
void extended_marching_tets(const Eigen::MatrixXd& tet_vertices,
                            const Eigen::MatrixXi& tet_connectivity,
                            const Eigen::VectorXd& isovalues,
                            double surface_isovalue,
                            Eigen::MatrixXd& tri_vertices,
                            Eigen::MatrixXi& tri_connectivity,
                            Eigen::VectorXi& tri_vertex_to_tet_vertex_map,
                            Eigen::MatrixXi& tri_vertex_to_tet_split_edge_map);

/// @brief Helper struct for alternative extended_marching_tets.
struct MarchingTetState
{
    /// List of triangle faces created (connectivity, origin tetrahedra)
    std::vector<std::pair<Eigen::RowVector3i, int>> faces;
    /// Set of unique faces (used to avoid creating duplicates)
    std::set<std::array<int, 3>> unique_faces;
    /// Vertices created at split edges
    std::vector<Eigen::RowVector3d> split_edge_vertices;
    /// Split edges for each corresponding entry in split_edge_vertices
    std::vector<std::array<int, 2>> split_edges;
    /// Note: the index will be negative to differentiate from existing vertices and new vertices
    /// from split edges. Corresponding vertex in split_edge_vertices is given by -1 * (value) - 1
    std::unordered_map<int64_t, int> split_edge_to_split_edge_vertex_map;
    /// Indicates if each of the input tetrahedra vertices should be included for the surface
    std::vector<bool> is_existing_vertex_included;
};

/// @brief Alternative marching tetrahedra algorithm that can be called multiple times before the
/// surface is recovered, allowing multiple surfaces to be processed creating a single result.
///
/// Notes:
/// - Each edge will only cut once, so if an edge is cut by call to this function, subsequent calls
/// with the same state will reuse the same point.
void extended_marching_tets(const Eigen::MatrixXd& tet_vertices,
                            const Eigen::MatrixXi& tet_connectivity,
                            const Eigen::VectorXd& isovalues,
                            const double& surface_isovalue,
                            MarchingTetState& state);

/// @brief Recovers isosurface from a call to the alternative extended_marching_tets.
void recover_isosurface(const Eigen::MatrixXd& tet_vertices,
                        const MarchingTetState& state,
                        Eigen::MatrixXd& tri_vertices,
                        Eigen::MatrixXi& tri_connectivity,
                        Eigen::VectorXi& tri_vertex_to_tet_vertex_map,
                        Eigen::MatrixXi& tri_vertex_to_tet_split_edge_map);


// New implementation to replace the old ones once complete

/// Implementation of marching tetrahedra that has two features beyond the standard implementation:
/// 1. A subset of the vertices may have an isovalue while the others do not, which is useful for
///    local operations like the medial surface method.
/// 2. Significantly more output information is provided and more cases were added to lookup tables
///    to track how the isosurface stiches into the larger tetrahedralization, which can be used
///    to remesh the tetrahedralization to accommodate the isosurface.
///
/// Note: there is an important consideration.  If a tetrahedra has some vertices with an isovalue
/// and and some vertices without an isovalue, then the isovalue shoule be 0 at each of those points
/// with an isovalue if you do not want the isosurface to terminate on a triangle face that is not
/// split during remeshing.  Of course, if you just want the isosurface and do not want to remesh
/// the tetredralization to accommodate it, then this restriction does not apply.
///
/// @param t - the tetrahedralization used to calculate the isovalues
/// @param vertices_with_isovalues - the list of vertices in `t` that have an isovalue, which can be
///        the entire set of vertices or a subset (see note above about an important consideration)
/// @param isovalues - the values from which to determine the isosurface corresponding to the list
///        of vertices `vertices_with_isovalues`
/// @param surface_isovalue - the value at which to find the surface (typically 0 for signed
///        distances)
///
/// @out surface_triangles - connectivity of the triangles for the oriented isosurface.  Indices
///      larger than those already in `t` indicate a vertex on an intersected edge with the ID given
///      by intersected_edge_ids[vertex id - t.vertices.size()]
/// @out intersected_edge_ids - the list of edges that were intersected by the isosurface (see
///      description of surface_triangles for more information)
/// @out point_on_split_edges - the coordinate of the point along the intersected edges
void extended_marching_tets(const Tetrahedralization& t,
                            const std::vector<mesh::LocalIndex>& vertices_with_isovalues,
                            const Eigen::VectorXd& isovalues,
                            const double surface_isovalue,
                            std::vector<std::array<mesh::LocalIndex, 3>>& surface_triangles,
                            std::vector<mesh::LocalIndex>& intersected_edge_ids,
                            std::vector<std::array<double, 3>>& point_on_intersected_edges);

/// @brief Helper struct for alternative extended_marching_tets.
struct MarchingTetState2
{
    /// List of triangle faces created (connectivity, origin tetrahedra)
    std::vector<std::pair<std::array<mesh::LocalIndex, 3>, int>> faces;
    /// Set of unique faces, which is used to quickly determine if a given face already exists (thus
    /// avoiding duplicate faces with successive calls).
    std::set<std::array<mesh::LocalIndex, 3>> unique_faces;
    /// Vertices created at split edges
    std::vector<std::array<double, 3>> split_edge_vertices;
    /// Length of all edges (uninitialized values take mesh::NOT_SET_LocalIndex)
    /// Indicates where in split_edge_vertices to look for the vertex along the split edge.
    std::vector<mesh::LocalIndex> edge_to_split_edge_vertex_map;
};

/// @out state - data that defines the isosurface(s) from successive calls to this function
/// @returns The list of indices for faces that define the surface from this call.
/*
std::vector<mesh::LocalIndex> extended_marching_tets(
  const Tetrahedralization& t,
  const std::vector<mesh::LocalIndex>& vertices_with_isovalues,
  const Eigen::VectorXd& isovalues,
  const double surface_isovalue,
  MarchingTetState2& state);
  */

void extended_marching_tets(const Tetrahedralization& t,
                            const std::vector<mesh::LocalIndex>& vertices_with_isovalues,
                            const Eigen::VectorXd& isovalues,
                            const double surface_isovalue,
                            MarchingTetState2& state);

/// @brief Recovers isosurface from a call to the alternative extended_marching_tets.
void recover_isosurface(const Tetrahedralization& t,
                        const MarchingTetState2& state,
                        std::vector<std::array<mesh::LocalIndex, 3>>& surface_triangles,
                        std::vector<mesh::LocalIndex>& intersected_edge_ids,
                        std::vector<std::array<double, 3>>& point_on_intersected_edges);

}