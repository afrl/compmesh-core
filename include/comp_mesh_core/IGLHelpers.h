#pragma once

#include "Tetrahedralization.h"
#include "Triangulation.h"
#include "igl/AABB.h"
#include "mesh/Mesh.h"

#include <utility>
#include <vector>

namespace composite_mesher
{

/// @brief A data structure holding the list of vertices and triangles for a surface and a
/// corresponding hierarchical bounding-box tree.
/// The vertices and triangles are stored in Eigen matrices as preferred in the IGL library, and the
/// hierarchical bounding-box tree allows fast distance and projection operations.
struct IGLTriSurf
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    igl::AABB<Eigen::MatrixXd, 3> tree;
};

struct IGLTetrahedralization
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi tetrahedra;
};

IGLTriSurf mesh_to_igl(const mesh::Mesh& m);

IGLTriSurf mesh_to_igl(const mesh::Mesh& m, const mesh::ComponentIndex& component_to_extract);

IGLTriSurf mesh_to_igl(const mesh::Mesh& m,
                       const std::vector<mesh::ComponentIndex>& components_to_extract);

IGLTriSurf triagulation_to_igl(const Triangulation& t, const Marker& marker_to_extract);

IGLTriSurf tetrahedralization_to_igl(const Tetrahedralization& t,
                                     const std::vector<mesh::LocalIndex>& triangles_to_extract);

/// Simple return struct for igl_to_mesh, @see igl_to_mesh
struct IGLToMeshResult
{
    /// Resulting mesh from IGL with duplicate triangles removed
    mesh::Mesh mesh;
    /// A map from the output element in mesh to the corresponding list of original triangles, which
    /// could be multiple due to duplicates.
    std::vector<std::vector<Eigen::Index>> element_to_origin_tris_map;
};

/// @brief Converts the Eigen matrix-based reprsentation of a triangle surface mesh (used by IGL) to
/// a Mesh, removing duplicates if any exist.
IGLToMeshResult igl_to_mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& triangles);

/// @brief Converts the Eigen matrix-based reprsentation of a triangle surface mesh (used by IGL) to
/// a Mesh, removing duplicates if any exist.
IGLToMeshResult igl_to_mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixX3i& triangles);

/// Simple return struct for igl_to_mesh, @see igl_to_mesh
struct IGLToTriangulationResult
{
    /// Resulting mesh from IGL with duplicate triangles removed
    Triangulation t;
    /// A map from the output element in mesh to the corresponding list of original triangles, which
    /// could be multiple due to duplicates.
    std::vector<std::vector<Eigen::Index>> element_to_origin_tris_map;
};

/// @brief Converts the Eigen matrix-based reprsentation of a triangle surface mesh (used by IGL) to
/// a Triangulation, removing duplicates if any exist.
IGLToTriangulationResult igl_to_triangulation(const Eigen::MatrixXd& vertices,
                                              const Eigen::MatrixXi& triangles);

double igl_sq_distance(const IGLTriSurf& igl_surf, const mesh::Point& p);

Eigen::Vector3d igl_closest_point(const IGLTriSurf& igl_surf, const Eigen::Vector3d& p);

/// Creates a list of the facets of the input surfaces the format used by IGL (including a
/// hierarchical bounding-box tree for fast distance and projection operations).  This is used
/// during optimizations to ensure nodes are projected back to the facets they should lie on.
std::vector<IGLTriSurf> decompose_into_igl_facets(const mesh::Mesh& combined_surface_mesh);

/// Same as the other decompose_into_igl_facets but for a Triangulation
std::vector<IGLTriSurf> decompose_into_igl_facets(const Triangulation& combined_surface_mesh);

void igl_closest_point_to_facet(const Eigen::MatrixX3d& vertices,
                                const IGLTriSurf& igl_surf,
                                Eigen::MatrixX3d& vectors_to_surface,
                                Eigen::VectorXd& distances_to_surface);

IGLTetrahedralization mesh_to_igl_tet(const mesh::Mesh& m);


double signed_volume(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
void orient_closed_surface(const Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool positive);
void orient_volume(const Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}