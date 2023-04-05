#pragma once

#include "mesh/Defines.h"
#include "mesh/Mesh.h"
#include "utility/DynamicArray.h"

#include <array>
#include <vector>

namespace composite_mesher
{
/// @brief Stores all information describing a tetrahedralizaiton with several maps needed during
/// remeshing.
struct Tetrahedralization
{
    // Each dynamic array contains pairs in the format: { global index, localized index }
    using LocalizedMap = std::vector<
      utility::DynamicArray<std::pair<mesh::LocalIndex, mesh::ConnIndex>, mesh::ConnIndex>>;

    /// Coordinates of all the vertices
    std::vector<std::array<double, 3>> vertices;
    /// All tetrahedra connectivity with indices corresponding to positions in vertices
    std::vector<std::array<mesh::LocalIndex, 4>> tetrahedra;
    /// All triangles connectivity with indices corresponding to positions in vertices
    std::vector<std::array<mesh::LocalIndex, 3>> triangles;
    /// All edges connectivity with indices corresponding to positions in vertices
    std::vector<std::array<mesh::LocalIndex, 2>> edges;

    // Edge-triangle relationships

    /// For each edge, this contains a list of pairs in the form:
    /// { triangle index, position of edge in the triangle }
    LocalizedMap edge_to_tris_map;
    /// List of edges in each triangle, with the order implying the position of the edge within
    /// the triangle
    std::vector<std::array<mesh::LocalIndex, 3>> tri_to_edges_map;

    // Edge-tetrahedra relationships

    /// For each edge, this contains a list of pairs in the form:
    /// { tet index, position of edge in the tet }
    LocalizedMap edge_to_tets_map;
    /// List of edges in each tetrahedron, with the order implying the position of the edge within
    /// the tetrahedron
    std::vector<std::array<mesh::LocalIndex, 6>> tet_to_edges_map;

    // Triangle-tetrahedra relationships

    /// List of triangles in each tetrahedron, with the order implying the position of the triangles
    /// within the tetrahedron
    std::vector<std::array<mesh::LocalIndex, 4>> tet_to_tri_map;
    /// The tetrahedra that contain the triangle as a face.  If the triangle is on the boundary then
    /// the value of the second will be mesh::NOT_SET_LocalIndex.
    std::vector<std::array<mesh::LocalIndex, 2>> tri_to_tet_map;
};

mesh::Mesh extract_surface(const Tetrahedralization& t);
mesh::Mesh extract_volume(const Tetrahedralization& t);

////////////////////////////////////////////
// Functions that are helpful for debugging

bool check_tri_edge_relations(const Tetrahedralization& t);
bool check_tet_edge_relations(const Tetrahedralization& t);
bool check_tet_tri_relations(const Tetrahedralization& t);

}