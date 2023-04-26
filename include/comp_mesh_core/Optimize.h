#pragma once

#include "mesh/Mesh.h"
#include "Triangulation.h"

namespace composite_mesher
{
/// @brief Optimizes a surface mesh using tangential smoothing, while preserving features.
void optimize_tangential_smoothing(mesh::Mesh& surface_mesh,
                                   const std::vector<std::array<mesh::LocalIndex, 2>>& edges,
                                   const mesh::LocalIndex& n_seams,
                                   const std::vector<mesh::LocalIndex>& edges_on_seams,
                                   const std::vector<mesh::LocalIndex>& edge_id_to_seam_id_map,
                                   const bool allow_smoothing_on_seams = true);

void optimize_taubin_smoothing(mesh::Mesh& surface_mesh,
                               const mesh::ComponentIndex& component_to_smooth,
                               const std::vector<bool>& is_node_constrained,
                               const double lambda              = 0.330,
                               const double mu                  = -0.331,
                               const int n_smoothing_iterations = 10);

void optimize_edge_flips(Triangulation& t, const int n_passes);
/*
// Approach #1

/// @brief Optimizes a surface mesh by refining the mesh, while preserving features.
void optimize_refinement(SurfaceMeshWithFeatures& feature_surface,
                         const double& target_edge_length);

void optimize_corefinement(SurfaceMeshWithFeatures& feature_surface,
                           const double target_mean_edge_length);


mesh::Mesh tetrahedralize(SurfaceMeshWithFeatures& feature_surface);

// Approach #2

/// Note: any facet_id for which is_facet_planar is true will be reduced to convex polygons (with
/// some of the polygons capturing holes)
///
/// Note: this doesn't really work.  If we let tetgen refine the surfaces, it creates far too much
/// refinement near sharp features. Additionally, it creates odd surface triangles that are very
/// sharp
mesh::Mesh tetrahedralize_via_direct_tetgen(const std::vector<bool>& is_facet_planar,
                                            SurfaceMeshWithFeatures& feature_surface);

// Approach #3

struct PLCFacet
{
    std::vector<std::vector<mesh::LocalIndex>> boundary_loops;
    std::vector<bool> is_boundary_loop_a_hole;
};

struct PLC
{
    std::vector<PLCFacet> facets;
};

PLC feature_surface_to_plc(const SurfaceMeshWithFeatures& feature_surface,
                           const std::vector<bool>& is_facet_planar);

void plc_to_polygonal_surface(const PLC& plc,
                              const SurfaceMeshWithFeatures& feature_surface,
                              const double max_triangle_area);

*/
}