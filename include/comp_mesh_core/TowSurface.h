#pragma once

#include "mesh/Mesh.h"

namespace composite_mesher
{
/// @brief A struct to store information about a tow
struct VTMSTowSurface
{
    enum class MeshType
    {
        POLY,
        TRI,
        TET
    };

    enum class NormalDirection
    {
        INWARD,
        OUTWARD
    };

    enum class SurfaceElementType
    {
        INTERIOR, ///< Polygon or triangle that is interior to the textile model (not clipped)
        BOUNDARY  ///< Polygon or triangle that lies on the boundary of the textile model (clipped)
    };

    int material_number = 0;
    std::vector<mesh::Point> medials;
    mesh::Mesh mesh;
    std::vector<SurfaceElementType> surface_element_types;
    NormalDirection surface_normal_direction = NormalDirection::INWARD;
    MeshType type                            = MeshType::POLY;
};

/// @brief Reads a .ctw file from VTMS into a mesh as a collection of polygonal elements that define
/// the bounds of a tow.
/// @param file_path - relative (to cwd) or absolute path to .ctw file
VTMSTowSurface read_vtms_ctw(const std::string& file_path);

/// @brief TODO document this
/// Note: does not read medials
VTMSTowSurface read_vtms_stw(const std::string& file_path);

/// @brief Creates a tow surface mesh from data structures instead of a .ctw file. @see
/// read_vtms_ctw
VTMSTowSurface create_tow_surface(const int material_number,
                               const std::vector<double>& points,
                               const std::vector<mesh::LocalIndex>& connectivities,
                               const std::vector<mesh::LocalIndex>& n_nodes_per_element,
                               const std::vector<double>& medials);

/// @brief Takes a polygonal discretization of a tow's surface and transforms it into a
///        quality triangulation of the tow's surface.
/// @param tow_surface - polygonal surface mesh for a tow returned by read_vtms_ctw
/// @returns Mean segment length
double triangulate_tow_surface(VTMSTowSurface& tow_surface, const bool refine = true);

/// @brief Write a triangulated tow surface to a .ply file, which can be read in by TetWild.
void write_to_ply(const mesh::Mesh& triangular_surface, const std::string& file_path);

}