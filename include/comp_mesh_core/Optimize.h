#pragma once

#include "Seams.h"
#include "Triangulation.h"
#include "mesh/Mesh.h"

namespace composite_mesher
{

/// @brief Optimizes a surface mesh using tangential smoothing, while preserving features.
// TODO delete
void optimize_tangential_smoothing(mesh::Mesh& surface_mesh,
                                   const std::vector<std::array<mesh::LocalIndex, 2>>& edges,
                                   const mesh::LocalIndex& n_seams,
                                   const std::vector<mesh::LocalIndex>& edges_on_seams,
                                   const std::vector<mesh::LocalIndex>& edge_id_to_seam_id_map,
                                   const bool allow_smoothing_on_seams = true);

/// @brief Optimizes a surface mesh using tangential smoothing, while preserving features.
void optimize_tangential_smoothing(Triangulation& t,
                                   const SeamInfo& seams,
                                   const int n_iterations,
                                   const bool allow_smoothing_on_seams = true);

void optimize_taubin_smoothing(mesh::Mesh& surface_mesh,
                               const mesh::ComponentIndex& component_to_smooth,
                               const std::vector<bool>& is_node_constrained,
                               const double lambda              = 0.330,
                               const double mu                  = -0.331,
                               const int n_smoothing_iterations = 10);

void optimize_edge_flips(Triangulation& t, const int n_passes);

}