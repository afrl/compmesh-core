#include "comp_mesh_core/VertexMap.h"

namespace composite_mesher
{
StaticVertexMapAccessor StaticVertexMap::get_indices(const Index& vertex_id) const
{
    return { *this, vertex_id };
}

StaticVertexMap build_vert_to_tri_map(const Triangulation& t)
{
    // First count the number of attached triangles for each vertex
    std::vector<int> n_tris_per_vertex(t.vertices.size(), 0);
    for (auto tri_id : utility::indices(t.triangles))
    {
        for (auto v_id : t.triangles[tri_id])
        {
            n_tris_per_vertex[v_id]++;
        }
    }

    StaticVertexMap m;
    m.offsets.resize(t.vertices.size() + 1, 0);
    for (auto v_id : utility::indices(t.vertices))
    {
        m.offsets[v_id + 1] = m.offsets[v_id] + n_tris_per_vertex[v_id];
    }

    m.id_list.resize(m.offsets.back());
    for (auto tri_id : utility::indices(t.triangles))
    {
        for (auto v_id : t.triangles[tri_id])
        {
            const auto position = m.offsets[v_id + 1] - n_tris_per_vertex[v_id];
            m.id_list[position] = tri_id;
            n_tris_per_vertex[v_id]--;
        }
    }

    return m;
}

}