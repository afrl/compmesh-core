#include "comp_mesh_core/CollapseEdges.h"

#include "simple_math/Utility.h"

namespace composite_mesher
{
namespace detail
{
struct CollapseState
{
    std::vector<Index> old_to_new_vertex_map;
    std::vector<Index> old_to_new_edge_map;
    std::vector<Index> old_to_new_triangle_map;
};

void collapse_edge(const Index& old_edge_id_to_collapse, Triangulation& t, CollapseState& state)
{
    using namespace simple_math::point_operators;

    auto new_edge_id_to_collapse = state.old_to_new_edge_map[old_edge_id_to_collapse];

    auto& edge_to_collapse = t.edges[new_edge_id_to_collapse];
    auto v1_id             = edge_to_collapse[0];
    auto v2_id             = edge_to_collapse[1];
    if (v1_id > v2_id)
    {
        std::swap(v1_id, v2_id);
    }
    t.vertices[v1_id] = 0.5 * (t.vertices[v1_id] + t.vertices[v2_id]);

    auto connected_tris = t.edge_to_tris_map[old_edge_id_to_collapse];
    for (auto tri : connected_tris)
    {
        // TODO delete tri
        auto tri_id = tri.first;
        auto edge_position_in_tri = tri.second;
        t.tri_to_edges_map[tri_id];
        t.triangles[tri_id];
    }

    state.old_to_new_edge_map[old_edge_id_to_collapse] = mesh::NOT_SET_LocalIndex;
}
}

CollapseResult collapse_edges(const std::vector<Index>& edges_to_collapse, Triangulation& t)
{
    detail::CollapseState state;
    // Initialize maps as 1-to-1 (unit map)
    state.old_to_new_vertex_map.resize(t.vertices.size());
    std::iota(state.old_to_new_vertex_map.begin(), state.old_to_new_vertex_map.end(), 0);
    state.old_to_new_edge_map.resize(t.edges.size());
    std::iota(state.old_to_new_edge_map.begin(), state.old_to_new_edge_map.end(), 0);
    state.old_to_new_triangle_map.resize(t.triangles.size());
    std::iota(state.old_to_new_triangle_map.begin(), state.old_to_new_triangle_map.end(), 0);

    throw _UTIL_MSG_EXCEPT("Not implemented yet.");
}

}