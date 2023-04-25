#pragma once

#include "Common.h"

#include "Triangulation.h"

namespace composite_mesher
{
/// @see segment_seams
struct SeamInfo
{
    Marker n_seams;
    // The list of edges that lie along a seam
    std::vector<Index> edges_on_seams;
    // A map from the edge ID to the seam ID.  Edges that do not lie on a seam will have the value
    // NOT_SET_SeamID.
    std::vector<Marker> edge_id_to_seam_id_map;
};


/// @brief Segemnts the seams that connect groups of triangles with different markers.
///
/// A seam refers to the edges that form the intersection of a set of surface patches with
/// differnet markers.  Segmentation is based on the combination of markers that intersect along the
/// seam.
SeamInfo segment_seams(const std::vector<std::array<Index, 2>>& edges,
                       const std::vector<std::array<Index, 3>>& tri_to_edges_map,
                       const Marker& n_markers,
                       const std::vector<Marker>& triangle_to_marker_map);

/// @brief Same as other segment_seams, but is convenient for Triangulation objects
SeamInfo segment_seams(const Triangulation& t);
}