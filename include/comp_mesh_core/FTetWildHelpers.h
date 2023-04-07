#pragma once

#include "mesh/Mesh.h"

namespace composite_mesher
{

/// @brief TODO
mesh::Mesh tetrahedralize_ftetwild(const mesh::Mesh& object_surface,
                                   const double ideal_relative_edge_length);

/// @brief TODO
mesh::Mesh tetrahedralize_ftetwild(const std::vector<mesh::Mesh>& objects_surfaces,
                                   const double ideal_relative_edge_length);


}