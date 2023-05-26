#pragma once

#include "comp_mesh_core/Common.h"
#include "mesh/Mesh.h"
#include "utility/DynamicArray.h"

#include <array>
#include <vector>

namespace composite_mesher
{
/// @brief Stores all information describing a triangulation of surfaces with several maps needed
/// during remeshing.
struct Triangulation
{
    // Each dynamic array contains pairs in the format: { global index, localized index }
    using LocalizedMap = std::vector<utility::DynamicArray<std::pair<Index, ConnIndex>, ConnIndex>>;

    /// Coordinates of all the vertices
    std::vector<std::array<double, 3>> vertices;
    /// All triangles connectivity with indices corresponding to positions in vertices
    std::vector<std::array<Index, 3>> triangles;
    /// Total number of markers, which should be max(tri_to_marker_map) + 1
    Marker n_markers = 0;
    /// A marker for each triangle
    std::vector<Marker> tri_to_marker_map;

    // Edges and edge-triangle relationships

    /// All edges connectivity with indices corresponding to positions in vertices
    std::vector<std::array<Index, 2>> edges;
    /// For each edge, this contains a list of pairs in the form:
    /// { triangle index, position of edge in the triangle }
    LocalizedMap edge_to_tris_map;
    /// List of edges in each triangle, with the order implying the position of the edge within
    /// the triangle
    std::vector<std::array<Index, 3>> tri_to_edges_map;
};

/// @brief Generates members of a Triangulation object that can be derived from others, including
/// - n_markers
/// - edges
/// - edge_to_tris_map
/// - tri_to_edges_map
void update_edges_and_maps(Triangulation& t);

struct EdgeLengthStats
{
    double min;
    double max;
    double mean;
    double variance;
    double std_dev;
};

EdgeLengthStats query_edge_length_stats(const Triangulation& t);

mesh::Mesh to_mesh(const Triangulation& t);
Triangulation from_mesh(const mesh::Mesh& m);
}