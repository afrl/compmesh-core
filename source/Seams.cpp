#include "comp_mesh_core/Seams.h"

#include "utility/DynamicBitset.h"

namespace composite_mesher
{
SeamInfo segment_seams(const std::vector<std::array<Index, 2>>& edges,
                       const std::vector<std::array<Index, 3>>& tri_to_edges_map,
                       const Marker& n_markers,
                       const std::vector<Marker>& triangle_to_marker_map)
{
    std::vector<utility::DynamicBitset<64>> edge_id_to_interface_bitset_map(
      edges.size(),
      utility::DynamicBitset<64>(n_markers));
    for (auto tri_id : utility::range(tri_to_edges_map.size()))
    {
        for (auto& edge_id : tri_to_edges_map[tri_id])
        {
            auto& tri_region_boundary = triangle_to_marker_map[tri_id];
            edge_id_to_interface_bitset_map[edge_id].set(tri_region_boundary, true);
        }
    }

    std::vector<Index> edges_on_seams;
    std::vector<Marker> edge_id_to_seam_id_map(edges.size(), NOT_SET_SeamID);
    std::map<utility::DynamicBitset<64>, std::size_t> bitset_to_seam_id_map;
    for (auto edge_id : utility::indices_as<Index>(edge_id_to_interface_bitset_map))
    {
        auto bitset = edge_id_to_interface_bitset_map[edge_id];
        if (bitset.count() > 1)
        {
            edges_on_seams.push_back(edge_id);

            auto find_iter = bitset_to_seam_id_map.find(bitset);
            if (find_iter == bitset_to_seam_id_map.end())
            {
                auto next_seam_id               = static_cast<Index>(bitset_to_seam_id_map.size());
                bitset_to_seam_id_map[bitset]   = next_seam_id;
                edge_id_to_seam_id_map[edge_id] = next_seam_id;
            }
            else
            {
                edge_id_to_seam_id_map[edge_id] = find_iter->second;
            }
        }
    }

    SeamInfo seam_info = { static_cast<Index>(bitset_to_seam_id_map.size()),
                           edges_on_seams,
                           edge_id_to_seam_id_map };

    return seam_info;
}

SeamInfo segment_seams(const Triangulation& t)
{
    return segment_seams(t.edges, t.tri_to_edges_map, t.n_markers, t.tri_to_marker_map);
}
}