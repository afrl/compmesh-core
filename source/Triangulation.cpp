#include "comp_mesh_core/Triangulation.h"

#include "comp_mesh_core/Common.h"
#include "mesh/algorithms/IO.h"
#include "mesh/common/ElementTypes.h"
#include "mesh/data/MeshAssociatedData.h"
#include "simple_math/Utility.h"

namespace composite_mesher
{
/// - n_markers
/// - edges
/// - edge_to_tris_map
/// - tri_to_edges_map
void update_edges_and_maps(Triangulation& t)
{
    t.n_markers = *std::max_element(t.tri_to_marker_map.begin(), t.tri_to_marker_map.end()) + 1;

    auto& edge_info = mesh::get_edge_info(mesh::ElementType::T2D3N);
    assert(edge_info.size() == 3);

    t.edges.clear();
    t.edge_to_tris_map.clear();
    t.tri_to_edges_map.resize(t.triangles.size());

    std::map<std::size_t, Index> edge_hash_map;
    for (auto tri_id : utility::indices(t.triangles))
    {
        for (auto i : utility::range(3))
        {
            assert(edge_info[i].type == mesh::EdgeType::E2N);

            auto tri_conn = t.triangles[tri_id];
            auto& edge_id = t.tri_to_edges_map[tri_id][i];

            auto& v1 = tri_conn[edge_info[i].connectivity[0]];
            auto& v2 = tri_conn[edge_info[i].connectivity[1]];

            auto hash      = hash_linear_edge(v1, v2);
            auto find_iter = edge_hash_map.find(hash);
            // This is a new edge
            if (find_iter == edge_hash_map.end())
            {
                auto edge_id = t.edges.size();
                t.edges.push_back({ v1, v2 });
                edge_hash_map[hash] = edge_id;

                t.edge_to_tris_map.push_back({ { tri_id, i } });
                t.tri_to_edges_map[tri_id][i] = edge_id;
            }
            // Edge has already been added
            else
            {
                auto edge_id = find_iter->second;
                t.edge_to_tris_map[edge_id].push_back({ tri_id, i });
                t.tri_to_edges_map[tri_id][i] = edge_id;
            }
        }
    }
}

EdgeLengthStats query_edge_length_stats(const Triangulation& t)
{
    // See https://www.johndcook.com/blog/standard_deviation/
    double mean     = 0.0;
    double variance = 0.0;
    double min      = std::numeric_limits<double>::max();
    double max      = std::numeric_limits<double>::min();
    for (auto k : utility::indices(t.edges))
    {
        auto& edge     = t.edges[k];
        auto length    = simple_math::distance(t.vertices[edge[0]], t.vertices[edge[1]]);
        auto prev_mean = mean;
        mean += (length - mean) / (double)(k + 1);
        if (k > 0)
        {
            variance += (length - prev_mean) * (length - mean) / (double)(k);
        }
        if (length < min)
        {
            min = length;
        }
        if (length > max)
        {
            max = length;
        }
    }
    return { min, max, mean, variance, sqrt(variance) };
}

mesh::Mesh to_mesh(const Triangulation& t)
{
    mesh::Mesh m;
    for (auto& p : t.vertices)
    {
        m.nodes().add(p);
    }
    for (auto i : utility::indices(t.triangles))
    {
        m.elements().add(mesh::ElementType::T2D3N,
                         { t.triangles[i][0], t.triangles[i][1], t.triangles[i][2] });
        m.elements().back().component_id() = t.tri_to_marker_map[i];
    }
    return m;
}
}