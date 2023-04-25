#include "comp_mesh_core/SplitEdges.h"

#include "comp_mesh_core/GlobalOptions.h"
#include "mesh/algorithms/IO.h"
#include "simple_math/Utility.h"

namespace composite_mesher
{

/// @brief Stores intermediate data for tetrahedralization routines.
struct UpdateState
{
    std::vector<mesh::LocalIndex> orig_edge_to_update_index;
    std::vector<mesh::LocalIndex> orig_tri_to_update_index;
    std::vector<mesh::LocalIndex> orig_tet_to_update_index;

    std::vector<mesh::LocalIndex> current_to_orig_edge_map;
    std::vector<mesh::LocalIndex> current_to_orig_tri_map;
    std::vector<mesh::LocalIndex> current_to_orig_tet_map;

    UpdateResult result;
};

namespace detail
{
struct TriangleNewEdgeInfo
{
    mesh::LocalIndex new_triangle_id;
    mesh::LocalIndex new_edge_id;
    bool is_edge_aligned;
};

/// @returns Map from original triangle -> { new triangle, new internal edge }
/// Note: the original triangle will always contain the first node of the original split edge,
/// and the triangle will contain the second node of the original split edge.
template<typename T>
std::map<mesh::LocalIndex, TriangleNewEdgeInfo> split_triangles_on_edge(
  const mesh::LocalIndex& orig_edge_id,
  const mesh::LocalIndex& new_split_edge_id,
  T& t)
{
    const auto& edges_info = mesh::get_edge_info(mesh::ElementType::T2D3N);

    auto& triangles        = t.triangles;
    auto& edges            = t.edges;
    auto& edge_to_tris_map = t.edge_to_tris_map;
    auto& tri_to_edges_map = t.tri_to_edges_map;

    auto n_1_id         = edges[orig_edge_id][0];
    auto new_node_index = edges[orig_edge_id][1];
    auto n_2_id         = edges[new_split_edge_id][1];

    // Add an entry for the new edge
    edge_to_tris_map.push_back({});

    // Map from original triangle -> { new triangle, new internal edge, is_edge_aligned }
    std::map<mesh::LocalIndex, TriangleNewEdgeInfo> tri_updates;

    // Note: parent_infos is intentionally not a reference since edge_id_to_triangles_ is modified
    // in the loop, which would invalidate an iterator
    auto parent_infos = edge_to_tris_map[orig_edge_id];
    for (auto& parent_info : parent_infos)
    {
        auto& parent_triangle_index    = parent_info.first;
        auto& split_edge_position_in_t = parent_info.second;
        auto& split_edge_info          = edges_info[split_edge_position_in_t];

        auto t_id     = parent_triangle_index;
        auto new_t_id = static_cast<mesh::LocalIndex>(triangles.size());
        triangles.push_back(triangles[t_id]);
        auto& t_conn     = triangles[t_id];
        auto& new_t_conn = triangles.back();

        // Add entries in the structures for the new triangle (the value of triangle_to_edge_map_
        // will be updated below)
        tri_to_edges_map.push_back(tri_to_edges_map[t_id]);

        // Find the third node (external_v) in the triangle (t) that is not on the edge being split
        mesh::ConnIndex external_conn_index = 0;
        for (auto i : utility::range(3))
        {
            if (i != split_edge_info.connectivity[0] && i != split_edge_info.connectivity[1])
            {
                external_conn_index = i;
            }
        }
        auto external_v = t_conn[external_conn_index];

        // Update connectivities
        t_conn[split_edge_info.connectivity[1]]     = new_node_index;
        new_t_conn[split_edge_info.connectivity[0]] = new_node_index;

        // So...this is confusing.  I worked out the logic, but it is very obscure and I don't
        // know a cleaner way to express it without lots of redundant code.
        //
        // Assume we are splitting edge 0-1:
        //
        //   original t        new_t    t
        //       o 2                 o 2
        //      / \                 /|\ 
        //     /   \               / | \           <-- triangle 0
        //    /     \             /  |n \ 
        // 1 o-------o 0       1 o---o---o 0
        //
        // Index of edge 0-n in t
        mesh::ConnIndex split_half_edge_position_in_t = split_edge_position_in_t;
        // Index of edge n-2 in t
        mesh::ConnIndex internal_edge_position_in_t = (split_edge_position_in_t + 1) % 3;
        // index of edge 2-0 in t
        mesh::ConnIndex external_edge_position_in_t = (split_edge_position_in_t + 2) % 3;
        // Index of edge n-1 in new_t
        mesh::ConnIndex split_half_edge_position_in_new_t = split_edge_position_in_t;
        // Index of edge 1-2 in new_t
        mesh::ConnIndex external_edge_position_in_new_t = (split_edge_position_in_t + 1) % 3;
        // index of edge 2-n in new_t
        mesh::ConnIndex internal_edge_position_in_new_t = (split_edge_position_in_t + 2) % 3;

        // This next piece of code updates the triangle and edge relationships.  However, it is
        // important to remember that this function is splitting global edge n_1->n_2, but the
        // direction of the edge might be opposite in current triangle being adjusted.  Furthermore,
        // the logic is a little different for the two cases (directions are the same vs directions
        // are opposite).
        //
        // If the direction of edge in the triangle matches the direction of the global edge, then:
        //  - triangle t gets the recycled edge (n_1->new_node), so no update needed
        //  - triangle new_t gets the new edge (new_node->n_2)
        auto edge_aligned = t_conn[split_edge_info.connectivity[0]] == n_1_id;
        if (edge_aligned)
        {
            // Update relationships relating to the new triangle for new_n->n_2
            tri_to_edges_map[new_t_id][split_half_edge_position_in_new_t] = new_split_edge_id;
            edge_to_tris_map[new_split_edge_id].push_back(
              { new_t_id, split_half_edge_position_in_new_t });
            tri_to_edges_map[new_t_id][split_half_edge_position_in_new_t] = new_split_edge_id;

            for (auto& tri_for_edge : edge_to_tris_map[orig_edge_id])
            {
                if (tri_for_edge.first == t_id)
                {
                    tri_for_edge.second = split_half_edge_position_in_t;
                }
            }
        }
        // If the direction of edge in the triangle is opposite of the direction of the global edge,
        // then:
        //  - triangle new_t gets the recycled edge (n_1->new_node), so no update needed
        //  - triangle t gets the new edge (new_node->n_2)
        else
        {
            // Update relationships relating to 't' for new_n->n_2
            tri_to_edges_map[t_id][split_half_edge_position_in_t] = new_split_edge_id;
            edge_to_tris_map[new_split_edge_id].push_back({ t_id, split_half_edge_position_in_t });
            tri_to_edges_map[t_id][split_half_edge_position_in_t] = new_split_edge_id;


            for (auto& tri_for_edge : edge_to_tris_map[orig_edge_id])
            {
                if (tri_for_edge.first == t_id)
                {
                    tri_for_edge.first  = new_t_id;
                    tri_for_edge.second = split_half_edge_position_in_new_t;
                }
            }
        }

        // Add the new internal edge new_n->external_n (guaranteed not to be a seam, even if all
        // edges were seams)
        auto new_internal_edge_id = static_cast<mesh::LocalIndex>(edge_to_tris_map.size());
        // Update edge maps
        edges.push_back({ new_node_index, external_v });
        edge_to_tris_map.push_back(
          { { t_id, internal_edge_position_in_t }, { new_t_id, internal_edge_position_in_new_t } });
        // Update triangle-to-edge map
        tri_to_edges_map[t_id][internal_edge_position_in_t]         = new_internal_edge_id;
        tri_to_edges_map[new_t_id][internal_edge_position_in_new_t] = new_internal_edge_id;

        // Update the maps for the external edge in the new triangle (since the parent triangle
        // index has changed for that edge)
        const auto& external_edge_id = tri_to_edges_map[new_t_id][external_edge_position_in_new_t];
        for (auto& parent_info : edge_to_tris_map[external_edge_id])
        {
            if (parent_info.first == t_id)
            {
                parent_info.first  = new_t_id;
                parent_info.second = external_edge_position_in_new_t;
            }
        }

        tri_updates[t_id] = { new_t_id, new_internal_edge_id, edge_aligned };
    }
    assert(tri_to_edges_map.size() == triangles.size());
    return tri_updates;
}

/// @returns Map from original tetrahedron -> new tetrahedron
struct TetrahedraUpdateResult
{
    struct NewInternalEdgeInfo
    {
        mesh::LocalIndex new_edge_id;
        mesh::LocalIndex new_node_id;
    };

    std::map<mesh::LocalIndex, mesh::LocalIndex> original_to_new_tet_map;
    std::vector<NewInternalEdgeInfo> new_internal_edges;
};

TetrahedraUpdateResult update_tetrahedra_on_split_edge(
  const mesh::LocalIndex& orig_edge_id,
  const mesh::LocalIndex& new_split_edge_id,
  const std::map<mesh::LocalIndex, TriangleNewEdgeInfo>& tri_updates,
  Tetrahedralization& t)
{
    auto& triangles        = t.triangles;
    auto& tets             = t.tetrahedra;
    auto& edges            = t.edges;
    auto& edge_to_tris_map = t.edge_to_tris_map;
    auto& edge_to_tets_map = t.edge_to_tets_map;
    auto& tet_to_edges_map = t.tet_to_edges_map;
    auto& tet_to_tri_map   = t.tet_to_tri_map;
    auto& tri_to_tet_map   = t.tri_to_tet_map;

    TetrahedraUpdateResult r;

    /////////////////////////////////////////////////////////////
    // Start reference data and tables

    // Reference Tetrahedron:
    //
    //  Vertices/Edges   Vertices/Faces
    //        v3               v3
    //        o                o
    //       /|\              /|\ 
    //    ↗e2 | e5↖          / | \<- f1 (back face)
    //     / e4↑ \          /f3|f2\ 
    // v0 o---|e1→o v2  v0 o---|---o v2
    //     \  |  /          \  |  /
    //    ↘e0 | e3↗          \ | /<- f0 (bottom face)
    //       \|/              \|/
    //        o                o
    //        v1               v1
    //
    // Reference Triangle:
    //
    //        v2
    //        o
    //       / \  
    //     e2   e1
    //     /     \ 
    // v0 o--e0---o v1

    // Given a split edge position within the tetrahedron, this table returns the indices of the
    // tri faces that were split to accomodate the split edge
    static const mesh::ConnIndex split_face_table[6][2] = { { 0, 3 }, { 0, 1 }, { 1, 3 },
                                                            { 0, 2 }, { 2, 3 }, { 1, 2 } };

    // Given a split edge position within the tetrahedron, this table returns the pair of edges that
    // are replaced with corresponding split face's internal edge (see split_face_table) for the
    // recycled tetrahedron keeps from the original tetrahedron.
    static const mesh::ConnIndex recycled_tet_internal_edge_replace_table[6][2] = {
        { 3, 4 }, { 3, 5 }, { 5, 4 }, { 1, 5 }, { 5, 2 }, { 2, 4 }
    };

    // Given a split edge position within the tetrahedron, this table returns the pair of edges that
    // are replaced with corresponding split face's internal edge (see split_face_table) for the
    // new tetrahedron keeps from the original tetrahedron.
    static const mesh::ConnIndex new_tet_internal_edge_replace_table[6][2] = { { 1, 2 }, { 0, 2 },
                                                                               { 1, 0 }, { 0, 4 },
                                                                               { 3, 0 }, { 1, 3 } };

    // Given a split edge position within the tetrahedron, this table returns the pair of edges that
    // the recycled tetrahedron keeps from the original tetrahedron.
    static auto& recycled_tet_edges_to_keep_table = new_tet_internal_edge_replace_table;

    // Given a split edge position within the tetrahedron, this table returns the pair of edges that
    // the new tetrahedron keeps from the original tetrahedron.
    static auto& new_tet_edges_to_keep_table = recycled_tet_internal_edge_replace_table;

    // Given a split edge position within the tetrahedron, this table returns the index of the
    // tri face the recycled tetrahedron keeps from the original tetrahedron.
    static const mesh::ConnIndex recycled_tet_face_to_keep[6] = { 1, 3, 0, 3, 0, 0 };

    // Given a split edge position within the tetrahedron, this table returns the index of the
    // tri face the new tetrahedron keeps from the original tetrahedron.
    static const mesh::ConnIndex new_tet_face_to_keep[6] = { 2, 2, 2, 1, 1, 3 };

    // The list of edges from the original tetrahedron that is common to both tetrahedra after the
    // split.
    static const mesh::ConnIndex common_orig_edge_table[6] = { 5, 4, 3, 2, 1, 0 };

    const auto& edge_info = mesh::get_edge_info(mesh::ElementType::T3D4N);
    const auto& face_info = mesh::get_face_info(mesh::ElementType::T3D4N);

    // End of reference data and tables
    /////////////////////////////////////////////////////////////


    auto n_1_id   = edges[orig_edge_id][0];
    auto n_new_id = edges[orig_edge_id][1];
    auto n_2_id   = edges[new_split_edge_id][1];

    // Note: parent_infos is intentionally not a reference since the map is modified
    // in the loop, which would invalidate an iterator
    auto orig_parent_tets = edge_to_tets_map[orig_edge_id];
    // This will be updated in the loop
    edge_to_tets_map[orig_edge_id].clear();
    // For the new split edge
    edge_to_tets_map.push_back({});

    tri_to_tet_map.resize(triangles.size());
    for (auto& u : tri_updates)
    {
        // For the split face
        tri_to_tet_map[u.second.new_triangle_id] = tri_to_tet_map[u.first];
        // For the internal edge of the split face
        edge_to_tets_map.push_back({});
    }

    assert(tri_to_tet_map.size() == triangles.size());
    assert(edge_to_tets_map.size() == edges.size());

    // Tracks the new edges for faces that are split that are also not a tracked surface triangle
    std::map<std::array<mesh::LocalIndex, 3>, mesh::LocalIndex> untracked_face_to_new_edge_map;

    for (auto& parent_info : orig_parent_tets)
    {
        auto& orig_tet_index             = parent_info.first;
        auto& split_edge_position_in_tet = parent_info.second;

        // Duplicate the connected tet
        auto new_tet_index = static_cast<mesh::LocalIndex>(tets.size());
        tets.push_back(tets[orig_tet_index]);
        auto& recycled_tet_conn = tets[orig_tet_index];
        auto& new_tet_conn      = tets.back();
        auto orig_tet_conn      = tets[orig_tet_index];
        auto orig_tet_edges     = tet_to_edges_map[orig_tet_index];

        r.original_to_new_tet_map[orig_tet_index] = new_tet_index;

        // Test if the edge direction follows the direction within the table given by
        // mesh::get_edge_info
        auto& split_edge_c = edge_info[split_edge_position_in_tet].connectivity;
        // Make sure the edge to tet map is not corrupted (ensure edge is in the right place)
        assert(orig_tet_conn[split_edge_c[0]] == n_1_id ||
               orig_tet_conn[split_edge_c[0]] == n_2_id);
        assert(orig_tet_conn[split_edge_c[1]] == n_1_id ||
               orig_tet_conn[split_edge_c[1]] == n_2_id);

        bool edge_aligned = orig_tet_conn[split_edge_c[0]] == n_1_id;

        // Change the connectivity for the two tetrahedra (resulting from the split)
        recycled_tet_conn[split_edge_c[1]] = n_new_id;
        new_tet_conn[split_edge_c[0]]      = n_new_id;

        // Reset tri-tet and edge-tet relationships
        auto orig_tet_to_tris          = tet_to_tri_map[orig_tet_index];
        tet_to_tri_map[orig_tet_index] = { mesh::NOT_SET_LocalIndex,
                                           mesh::NOT_SET_LocalIndex,
                                           mesh::NOT_SET_LocalIndex,
                                           mesh::NOT_SET_LocalIndex };
        assert(tet_to_tri_map.size() == new_tet_index);
        tet_to_tri_map.push_back({ mesh::NOT_SET_LocalIndex,
                                   mesh::NOT_SET_LocalIndex,
                                   mesh::NOT_SET_LocalIndex,
                                   mesh::NOT_SET_LocalIndex });
        tet_to_edges_map[orig_tet_index] = { mesh::NOT_SET_LocalIndex, mesh::NOT_SET_LocalIndex,
                                             mesh::NOT_SET_LocalIndex, mesh::NOT_SET_LocalIndex,
                                             mesh::NOT_SET_LocalIndex, mesh::NOT_SET_LocalIndex };
        assert(tet_to_edges_map.size() == new_tet_index);
        tet_to_edges_map.push_back({ mesh::NOT_SET_LocalIndex,
                                     mesh::NOT_SET_LocalIndex,
                                     mesh::NOT_SET_LocalIndex,
                                     mesh::NOT_SET_LocalIndex,
                                     mesh::NOT_SET_LocalIndex,
                                     mesh::NOT_SET_LocalIndex });

        // Add a relation for the common edge to the new tetrahedron (it already has a relation to
        // the original tetrahedron)
        auto& common_edge_position = common_orig_edge_table[split_edge_position_in_tet];
        auto& common_edge          = orig_tet_edges[common_edge_position];
        tet_to_edges_map[orig_tet_index][common_edge_position] = common_edge;
        tet_to_edges_map[new_tet_index][common_edge_position]  = common_edge;
        edge_to_tets_map[common_edge].push_back(
          { new_tet_index, common_orig_edge_table[split_edge_position_in_tet] });

        // Update tet-edge relationships for the split edge
        if (edge_aligned)
        {
            tet_to_edges_map[orig_tet_index][split_edge_position_in_tet] = orig_edge_id;
            tet_to_edges_map[new_tet_index][split_edge_position_in_tet]  = new_split_edge_id;
            edge_to_tets_map[orig_edge_id].push_back(
              { orig_tet_index, split_edge_position_in_tet });
            edge_to_tets_map[new_split_edge_id].push_back(
              { new_tet_index, split_edge_position_in_tet });
        }
        else
        {
            tet_to_edges_map[new_tet_index][split_edge_position_in_tet]  = orig_edge_id;
            tet_to_edges_map[orig_tet_index][split_edge_position_in_tet] = new_split_edge_id;
            edge_to_tets_map[orig_edge_id].push_back({ new_tet_index, split_edge_position_in_tet });
            edge_to_tets_map[new_split_edge_id].push_back(
              { orig_tet_index, split_edge_position_in_tet });
        }

        // Update tet-tri relationships for the unsplit faces that are kept from the original
        // tetrahedron.
        tet_to_tri_map[orig_tet_index][recycled_tet_face_to_keep[split_edge_position_in_tet]] =
          orig_tet_to_tris[recycled_tet_face_to_keep[split_edge_position_in_tet]];
        // Note: the corresponding tri_to_tet_map already contains the original tet, and the
        // position did not change
        auto& face_for_new_tet = orig_tet_to_tris[new_tet_face_to_keep[split_edge_position_in_tet]];
        tet_to_tri_map[new_tet_index][new_tet_face_to_keep[split_edge_position_in_tet]] =
          face_for_new_tet;
        if (face_for_new_tet != mesh::NOT_SET_LocalIndex)
        {
            for (auto& t : tri_to_tet_map[face_for_new_tet])
            {
                if (t == orig_tet_index)
                {
                    t = new_tet_index;
                }
            }
        }

        // Update the tet-edge relationships for the edges that are kept from the original tet
        for (auto& edge_position : recycled_tet_edges_to_keep_table[split_edge_position_in_tet])
        {
            auto& edge_id = orig_tet_edges[edge_position];
            // edge_to_tets_map[edge_id] already has a reference to orig_tet_index
            tet_to_edges_map[orig_tet_index][edge_position] = edge_id;
        }
        for (auto& edge_position : new_tet_edges_to_keep_table[split_edge_position_in_tet])
        {
            auto& edge_id = orig_tet_edges[edge_position];
            for (auto& edge_to_tet_info : edge_to_tets_map[edge_id])
            {
                if (edge_to_tet_info.first == orig_tet_index)
                {
                    edge_to_tet_info.first  = new_tet_index;
                    edge_to_tet_info.second = edge_position;
                }
            }
            tet_to_edges_map[new_tet_index][edge_position] = edge_id;
        }
        // TODO check that push_back is correct rather than replacement

        // Address the split triangle faces and the new internal edge
        // Loop over the split faces (order does not matter)
        for (auto i : utility::range(2))
        {
            auto split_face_position = split_face_table[split_edge_position_in_tet][i];

            mesh::LocalIndex new_internal_edge_index = mesh::NOT_SET_LocalIndex;

            auto orig_tri_that_was_split = orig_tet_to_tris[split_face_position];
            if (orig_tri_that_was_split == mesh::NOT_SET_LocalIndex)
            {
                // Though the face is not a triangle that is being tracked, there still needs to be
                // a new internal edge along the face due to splitting the current edge.

                // Create hash of the face
                std::array<mesh::LocalIndex, 3> tri_face;
                auto& split_face_info = face_info[split_face_position];
                for (auto j : utility::range(3))
                {
                    tri_face[j] = orig_tet_conn[split_face_info.connectivity[j]];
                }
                std::sort(tri_face.begin(), tri_face.end());

                auto find_iter = untracked_face_to_new_edge_map.find(tri_face);
                if (find_iter == untracked_face_to_new_edge_map.end())
                {
                    // An internal edge has not been created for the untracked face yet
                    // Find the 3rd node on the triangle that should connect to the new node
                    mesh::LocalIndex third_node = mesh::NOT_SET_LocalIndex;
                    for (auto& v : tri_face)
                    {
                        if (v != n_1_id && v != n_2_id)
                        {
                            assert(third_node == mesh::NOT_SET_LocalIndex);
                            third_node = v;
                        }
                    }
                    assert(third_node != mesh::NOT_SET_LocalIndex);

                    new_internal_edge_index = edges.size();
                    edges.push_back({ third_node, n_new_id });
                    edge_to_tets_map.push_back({});
                    // There is no tracked triangles with this edge (since the program made it
                    // here), but it needs an entry
                    edge_to_tris_map.push_back({});
                    r.new_internal_edges.push_back({ new_internal_edge_index, n_new_id });
                    untracked_face_to_new_edge_map[tri_face] = new_internal_edge_index;
                }
                else
                {
                    new_internal_edge_index = find_iter->second;
                }
                assert(new_internal_edge_index != mesh::NOT_SET_LocalIndex);
            }
            else
            {
                const auto& orig_split_tri_index = orig_tri_that_was_split;
                auto find_iter                   = tri_updates.find(orig_split_tri_index);
                if (find_iter == tri_updates.end())
                {
                    throw _UTIL_MSG_EXCEPT(
                      "Something went wrong when updating the tetrahedra for the "
                      "split edge.  The lookup tables might be wrong.");
                }
                const auto& new_split_tri_index = find_iter->second.new_triangle_id;
                new_internal_edge_index         = find_iter->second.new_edge_id;
                const auto& is_tri_edge_aligned = find_iter->second.is_edge_aligned;

                // Update tet-tri relationships for the split triangles
                if (edge_aligned == is_tri_edge_aligned)
                {
                    // The recycled tet will get the original tri
                    // The new tet will get the new tri created from the split
                    if (tri_to_tet_map[orig_split_tri_index][0] == orig_tet_index)
                    {
                        assert(tri_to_tet_map[new_split_tri_index][0] == orig_tet_index);
                        tri_to_tet_map[new_split_tri_index][0] = new_tet_index;
                    }
                    else if (tri_to_tet_map[orig_split_tri_index][1] == orig_tet_index)
                    {
                        assert(tri_to_tet_map[new_split_tri_index][1] == orig_tet_index);
                        tri_to_tet_map[new_split_tri_index][1] = new_tet_index;
                    }
                    tet_to_tri_map[new_tet_index][split_face_position]  = new_split_tri_index;
                    tet_to_tri_map[orig_tet_index][split_face_position] = orig_split_tri_index;
                }
                else
                {
                    // The new tet will get the original tri
                    // The recycled tet will get the new tri created from the split
                    if (tri_to_tet_map[orig_split_tri_index][0] == orig_tet_index)
                    {
                        assert(tri_to_tet_map[new_split_tri_index][0] == orig_tet_index);
                        tri_to_tet_map[orig_split_tri_index][0] = new_tet_index;
                    }
                    else if (tri_to_tet_map[orig_split_tri_index][1] == orig_tet_index)
                    {
                        assert(tri_to_tet_map[new_split_tri_index][1] == orig_tet_index);
                        tri_to_tet_map[orig_split_tri_index][1] = new_tet_index;
                    }
                    tet_to_tri_map[new_tet_index][split_face_position]  = orig_split_tri_index;
                    tet_to_tri_map[orig_tet_index][split_face_position] = new_split_tri_index;
                }
            }

            // Update tet-edge relationships for the new internal edge
            auto& edge_to_replace_in_recycled_tet =
              recycled_tet_internal_edge_replace_table[split_edge_position_in_tet][i];
            tet_to_edges_map[orig_tet_index][edge_to_replace_in_recycled_tet] =
              new_internal_edge_index;
            edge_to_tets_map[new_internal_edge_index].push_back(
              { orig_tet_index, edge_to_replace_in_recycled_tet });

            // Update tet-edge relationships for the new internal edge
            auto& edge_to_replace_in_new_tet =
              new_tet_internal_edge_replace_table[split_edge_position_in_tet][i];
            tet_to_edges_map[new_tet_index][edge_to_replace_in_new_tet] = new_internal_edge_index;
            edge_to_tets_map[new_internal_edge_index].push_back(
              { new_tet_index, edge_to_replace_in_new_tet });
        }
    }

    return r;
}

void apply_tri_updates(const mesh::LocalIndex& new_vertex_id,
                       const std::map<mesh::LocalIndex, TriangleNewEdgeInfo>& tri_updates,
                       UpdateState& state)
{
    state.current_to_orig_edge_map.resize(state.current_to_orig_edge_map.size() +
                                            tri_updates.size(),
                                          mesh::NOT_SET_LocalIndex);

    for (auto update : tri_updates)
    {
        auto new_internal_edge_id = update.second.new_edge_id;
        state.result.updated_edges.push_back({ UpdatedEdgeInfo::NEW_EDGE,
                                               mesh::NOT_SET_LocalIndex,
                                               new_internal_edge_id,
                                               new_vertex_id });
    }

    state.current_to_orig_tri_map.resize(state.current_to_orig_tri_map.size() + tri_updates.size(),
                                         mesh::NOT_SET_LocalIndex);
    for (auto update : tri_updates)
    {
        if (state.current_to_orig_tri_map[update.first] != mesh::NOT_SET_LocalIndex)
        {
            // This triangle was the product of a previous update
            const auto orig_tri_id = state.current_to_orig_tri_map[update.first];
            const auto new_tri_id  = update.second.new_triangle_id;
            // If there is a current triangle created from an original triangle, then there has to
            // be an update entry
            assert(state.orig_tri_to_update_index[orig_tri_id] != mesh::NOT_SET_LocalIndex);
            // The original triangle already has an update from another edge split
            const auto update_index = state.orig_tri_to_update_index[orig_tri_id];
            state.result.updated_triangles[update_index].new_tri_ids.push_back(new_tri_id);
            state.current_to_orig_tri_map[new_tri_id] = orig_tri_id;
        }
        else
        {
            const auto orig_tri_id = update.first;
            const auto new_tri_id  = update.second.new_triangle_id;
            if (state.orig_tri_to_update_index[orig_tri_id] == mesh::NOT_SET_LocalIndex)
            {
                // There has not been an update entry for this triangle
                const auto update_index                     = state.result.updated_triangles.size();
                state.orig_tri_to_update_index[orig_tri_id] = update_index;
                state.current_to_orig_tri_map[new_tri_id]   = orig_tri_id;
                state.result.updated_triangles.push_back({ orig_tri_id, { new_tri_id } });
            }
            else
            {
                // There was already an update entry for this triangle, so just append to it
                const auto update_index = state.orig_tri_to_update_index[orig_tri_id];
                state.current_to_orig_tri_map[new_tri_id] = orig_tri_id;
                state.result.updated_triangles[update_index].new_tri_ids.push_back(new_tri_id);
            }
        }
    }
}

void apply_tet_updates(const TetrahedraUpdateResult& tet_updates, UpdateState& state)
{

    state.current_to_orig_edge_map.resize(state.current_to_orig_edge_map.size() +
                                            tet_updates.new_internal_edges.size(),
                                          mesh::NOT_SET_LocalIndex);

    for (auto update : tet_updates.new_internal_edges)
    {
        auto new_internal_edge_id = update.new_edge_id;
        auto new_node_id          = update.new_node_id;
        state.result.updated_edges.push_back({ UpdatedEdgeInfo::NEW_EDGE,
                                               mesh::NOT_SET_LocalIndex,
                                               new_internal_edge_id,
                                               new_node_id });
    }
    state.current_to_orig_tet_map.resize(state.current_to_orig_tet_map.size() +
                                           tet_updates.original_to_new_tet_map.size(),
                                         mesh::NOT_SET_LocalIndex);
    for (auto update : tet_updates.original_to_new_tet_map)
    {
        if (state.current_to_orig_tet_map[update.first] != mesh::NOT_SET_LocalIndex)
        {
            // This tetrahedron was the product of a previous update
            const auto orig_tet_id = state.current_to_orig_tet_map[update.first];
            const auto new_tet_id  = update.second;
            // If there is a current tetrahedron created from an original tetrahedron, then there
            // has to be an update entry
            assert(state.orig_tet_to_update_index[orig_tet_id] != mesh::NOT_SET_LocalIndex);
            // The original tetrahedron already has an update from another edge split
            const auto update_index = state.orig_tet_to_update_index[orig_tet_id];
            state.result.updated_tetrahedra[update_index].new_tet_ids.push_back(new_tet_id);
            state.current_to_orig_tet_map[new_tet_id] = orig_tet_id;
        }
        else
        {
            const auto orig_tet_id = update.first;
            const auto new_tet_id  = update.second;
            if (state.orig_tet_to_update_index[orig_tet_id] == mesh::NOT_SET_LocalIndex)
            {
                // There has not been an update entry for this tetrahedron
                const auto update_index = state.result.updated_tetrahedra.size();
                state.orig_tet_to_update_index[orig_tet_id] = update_index;
                state.current_to_orig_tet_map[new_tet_id]   = orig_tet_id;
                state.result.updated_tetrahedra.push_back({ orig_tet_id, { new_tet_id } });
            }
            else
            {
                // There was already an update entry for this tetrahedron, so just append to it
                const auto update_index = state.orig_tet_to_update_index[orig_tet_id];
                state.current_to_orig_tet_map[new_tet_id] = orig_tet_id;
                state.result.updated_tetrahedra[update_index].new_tet_ids.push_back(new_tet_id);
            }
        }
    }
}

/// @param edge_id - the index of the edge to split
/// @param existing_vertex_id_for_split - the index of the vertex that has already been created for
/// the split edge
void split_single_edge(const mesh::LocalIndex& edge_id,
                       const mesh::LocalIndex& existing_vertex_id_for_split,
                       Tetrahedralization& t,
                       UpdateState& state)
{
    using namespace simple_math::point_operators;

    // Not allowed to split the same edge twice
    assert(state.orig_edge_to_update_index[edge_id] == mesh::NOT_SET_LocalIndex);

    auto& vertices = t.vertices;
    auto& edges    = t.edges;

    auto v_1 = edges[edge_id][0];
    auto v_2 = edges[edge_id][1];
    auto p_1 = vertices[v_1];
    auto p_2 = vertices[v_2];

    mesh::LocalIndex new_v = existing_vertex_id_for_split;
    if (new_v == mesh::NOT_SET_LocalIndex)
    {
        new_v          = vertices.size();
        auto mid_point = 0.5 * (p_1 + p_2);
        vertices.push_back(mid_point);
    }


    // Create an ID for new_n->n_2, which is the new edge along the split edge
    // (n_1->new_n takes the place of the existing edge)
    auto new_n_n_2_edge_id = edges.size();
    edges.push_back({ new_v, v_2 });
    edges[edge_id][1] = new_v;

    auto tri_updates = detail::split_triangles_on_edge(edge_id, new_n_n_2_edge_id, t);

    if (GlobalOptions::get().check_common_errors)
    {
        assert(check_tri_edge_relations(t));
    }

    auto tet_updates =
      detail::update_tetrahedra_on_split_edge(edge_id, new_n_n_2_edge_id, tri_updates, t);

    if (GlobalOptions::get().check_common_errors)
    {
        assert(check_tet_edge_relations(t));
        assert(check_tri_edge_relations(t));
        assert(check_tet_tri_relations(t));
    }


    apply_tri_updates(new_v, tri_updates, state);
    apply_tet_updates(tet_updates, state);

    return;
}

void split_single_edge(const mesh::LocalIndex& edge_id,
                       const mesh::LocalIndex& existing_vertex_id_for_split,
                       Triangulation& t,
                       UpdateState& state)
{
    using namespace simple_math::point_operators;

    // Not allowed to split the same edge twice
    assert(state.orig_edge_to_update_index[edge_id] == mesh::NOT_SET_LocalIndex);

    auto& vertices = t.vertices;
    auto& edges    = t.edges;

    auto v_1 = edges[edge_id][0];
    auto v_2 = edges[edge_id][1];
    auto p_1 = vertices[v_1];
    auto p_2 = vertices[v_2];

    mesh::LocalIndex new_v = existing_vertex_id_for_split;
    if (new_v == mesh::NOT_SET_LocalIndex)
    {
        new_v          = vertices.size();
        auto mid_point = 0.5 * (p_1 + p_2);
        vertices.push_back(mid_point);
    }


    // Create an ID for new_n->n_2, which is the new edge along the split edge
    // (n_1->new_n takes the place of the existing edge)
    auto new_n_n_2_edge_id = edges.size();
    edges.push_back({ new_v, v_2 });
    edges[edge_id][1] = new_v;

    auto tri_updates = detail::split_triangles_on_edge(edge_id, new_n_n_2_edge_id, t);

    if (GlobalOptions::get().check_common_errors)
    {
        assert(check_tri_edge_relations(t));
    }

    apply_tri_updates(new_v, tri_updates, state);

    return;
}

void update_triangulation_markers(const UpdateResult& update_result, Triangulation& t)
{
    t.tri_to_marker_map.resize(t.triangles.size(), NOT_SET_Marker);
    for (auto& tri_update : update_result.updated_triangles)
    {
        auto& orig_marker = t.tri_to_marker_map[tri_update.orig_tri_id];
        for (auto& new_tri_id : tri_update.new_tri_ids)
        {
            t.tri_to_marker_map[new_tri_id] = orig_marker;
        }
    }
}

} // end namespace detail


UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split, Tetrahedralization& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tri_to_update_index.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tet_to_update_index.resize(t.tetrahedra.size(), mesh::NOT_SET_LocalIndex);

    state.current_to_orig_edge_map.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tri_map.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tet_map.resize(t.tetrahedra.size(), mesh::NOT_SET_LocalIndex);

    assert(t.tri_to_edges_map.size() == t.triangles.size());

    for (auto edge_id : edges_to_split)
    {
        detail::split_single_edge(edge_id, mesh::NOT_SET_LocalIndex, t, state);

        // auto m = extract_volume(t);
        // mesh::write_vtk(m, "volume_after_splitting_" + std::to_string(edge_id) + ".vtk");
    }

    return state.result;
}

UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         const std::vector<mesh::LocalIndex>& existing_vertices_for_split_edges,
                         Tetrahedralization& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tri_to_update_index.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tet_to_update_index.resize(t.tetrahedra.size(), mesh::NOT_SET_LocalIndex);

    state.current_to_orig_edge_map.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tri_map.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tet_map.resize(t.tetrahedra.size(), mesh::NOT_SET_LocalIndex);

    for (auto i : utility::indices(edges_to_split))
    {
        detail::split_single_edge(edges_to_split[i],
                                  existing_vertices_for_split_edges[i],
                                  t,
                                  state);

        // auto m = extract_volume(t);
        // mesh::write_vtk(m, "volume_after_splitting_" + std::to_string(edges_to_split[i]) +
        // ".vtk");
    }

    return state.result;
}

UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split, Triangulation& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tri_to_update_index.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);

    state.current_to_orig_edge_map.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tri_map.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);

    assert(t.tri_to_edges_map.size() == t.triangles.size());

    for (auto edge_id : edges_to_split)
    {
        detail::split_single_edge(edge_id, mesh::NOT_SET_LocalIndex, t, state);

        // auto m = extract_volume(t);
        // mesh::write_vtk(m, "volume_after_splitting_" + std::to_string(edge_id) + ".vtk");
    }

    detail::update_triangulation_markers(state.result, t);

    return state.result;
}

UpdateResult split_edges(const std::vector<mesh::LocalIndex>& edges_to_split,
                         const std::vector<mesh::LocalIndex>& existing_vertices_for_split_edges,
                         Triangulation& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.orig_tri_to_update_index.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);

    state.current_to_orig_edge_map.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    state.current_to_orig_tri_map.resize(t.triangles.size(), mesh::NOT_SET_LocalIndex);

    for (auto i : utility::indices(edges_to_split))
    {
        detail::split_single_edge(edges_to_split[i],
                                  existing_vertices_for_split_edges[i],
                                  t,
                                  state);

        // auto m = extract_volume(t);
        // mesh::write_vtk(m, "volume_after_splitting_" + std::to_string(edges_to_split[i]) +
        // ".vtk");
    }

    detail::update_triangulation_markers(state.result, t);

    return state.result;
}

UpdateResult split_constrained_edges(Triangulation& t, SeamInfo& seams)
{
    std::vector<bool> is_node_on_seam(t.vertices.size(), false);
    for (auto edge_id : seams.edges_on_seams)
    {
        auto n1_id             = t.edges[edge_id][0];
        is_node_on_seam[n1_id] = true;
        auto n2_id             = t.edges[edge_id][1];
        is_node_on_seam[n2_id] = true;
    }

    std::vector<mesh::LocalIndex> edges_to_split;
    for (auto edge_id : utility::indices(t.edges))
    {
        auto v_1                     = t.edges[edge_id][0];
        auto v_2                     = t.edges[edge_id][1];
        bool edge_points_constrained = is_node_on_seam[v_1] && is_node_on_seam[v_2];

        if (edge_points_constrained && seams.edge_id_to_seam_id_map[edge_id] == NOT_SET_SeamID)
        {
            // The edge is not a seam (it is a free edge), so it should be split
            edges_to_split.push_back(edge_id);
        }
    }
    std::cout << "* split_constrained_edges - Splitting " << edges_to_split.size() << " edges."
              << std::endl;

    auto update_result = split_edges(edges_to_split, t);

    // Split edges will update the Triangulation structure, but now update seam info
    seams.edge_id_to_seam_id_map.resize(t.edges.size(), NOT_SET_SeamID);
    for (auto& update : update_result.updated_edges)
    {
        if (update.update_type == UpdatedEdgeInfo::SUBDIVIDED_EDGE)
        {
            seams.edges_on_seams.push_back(update.new_edge_id);
            seams.edge_id_to_seam_id_map[update.new_edge_id] =
              seams.edge_id_to_seam_id_map[update.orig_edge_id];
        }
    }

    return update_result;
}
}