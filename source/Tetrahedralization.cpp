#include "comp_mesh_core/Tetrahedralization.h"

#include "comp_mesh_core/Common.h"
#include "mesh/algorithms/IO.h"
#include "mesh/common/ElementTypes.h"
#include "mesh/data/MeshAssociatedData.h"

namespace composite_mesher
{
mesh::Mesh extract_surface(const Tetrahedralization& t)
{
    mesh::Mesh m;
    for (auto& p : t.vertices)
    {
        m.nodes().add(p);
    }
    for (auto& t : t.triangles)
    {
        m.elements().add(mesh::ElementType::T2D3N, { t[0], t[1], t[2] });
    }
    return m;
}

mesh::Mesh extract_volume(const Tetrahedralization& t)
{
    mesh::Mesh m;
    for (auto& p : t.vertices)
    {
        m.nodes().add(p);
    }
    for (auto& t : t.tetrahedra)
    {
        m.elements().add(mesh::ElementType::T3D4N, { t[0], t[1], t[2], t[3] });
    }
    return m;
}

void write_debug_mesh_state_to_vtk(const Tetrahedralization& t, const std::string& filestem)
{
    {
        auto m = extract_surface(t);

        mesh::MeshAssociatedData data;

        if (!t.tri_to_tet_map.empty())
        {
            data.add<mesh::DataFieldType::ELEMENT>(
              mesh::DataField<mesh::LocalIndex>::cref("tri_to_tet_map",
                                                      { "tet_1", "tet_2" },
                                                      t.tri_to_tet_map[0].data(),
                                                      t.tri_to_tet_map.size(),
                                                      t.tri_to_tet_map[0].size()));
        }

        data.add<mesh::DataFieldType::ELEMENT>(
          mesh::DataField<mesh::LocalIndex>::cref("tri_to_edges_map",
                                                  { "edge_1", "edge_2", "edge_3" },
                                                  t.tri_to_edges_map[0].data(),
                                                  t.tri_to_edges_map.size(),
                                                  t.tri_to_edges_map[0].size()));

        mesh::write_vtk(m, filestem + "_triangles.vtk", data);
    }

    {
        auto m = extract_volume(t);

        mesh::MeshAssociatedData data;

        data.add<mesh::DataFieldType::ELEMENT>(
          mesh::DataField<mesh::LocalIndex>::cref("tet_to_tri_map",
                                                  { "tri_1", "tri_2", "tri_3", "tri_4" },
                                                  t.tet_to_tri_map[0].data(),
                                                  t.tet_to_tri_map.size(),
                                                  t.tet_to_tri_map[0].size()));

        data.add<mesh::DataFieldType::ELEMENT>(mesh::DataField<mesh::LocalIndex>::cref(
          "tet_to_edges_map",
          { "edge_1", "edge_2", "edge_3", "edge_4", "edge_5", "edge_6" },
          t.tet_to_edges_map[0].data(),
          t.tet_to_edges_map.size(),
          t.tet_to_edges_map[0].size()));

        mesh::write_vtk(m, filestem + "_tetrahera.vtk", data);
    }
}

template<typename T>
bool check_tri_edge_relations_kernel(const T& t)
{
    auto& triangles        = t.triangles;
    auto& edges            = t.edges;
    auto& vertices         = t.vertices;
    auto& tri_to_edges_map = t.tri_to_edges_map;
    auto& edge_to_tris_map = t.edge_to_tris_map;

    // Verify triangle-edge relationship maps
    for (auto tri_id : utility::indices(triangles))
    {
        for (auto i : utility::range(3))
        {
            auto& v_id = triangles[tri_id][i];
            if (!(v_id >= 0 && v_id < vertices.size()))
            {
                std::cout << "Triangle " + std::to_string(tri_id) + " references invalid vertex " +
                               std::to_string(v_id)
                          << std::endl;
            }
            _UTIL_ASSERT(v_id >= 0 && v_id < vertices.size(),
                         "Triangle " + std::to_string(tri_id) + " references invalid vertex " +
                           std::to_string(v_id));
        }

        auto& edge_info = mesh::get_edge_info(mesh::ElementType::T2D3N);
        assert(edge_info.size() == 3);

        // Loop over edges in the triangle
        for (auto i : utility::range(edge_info.size()))
        {
            assert(edge_info[i].type == mesh::EdgeType::E2N);

            auto tri_conn = triangles[tri_id];
            auto& edge_id = tri_to_edges_map[tri_id][i];

            // Check that every vertex in the edge exists in the triangle
            for (auto v : edges[edge_id])
            {
                _UTIL_ASSERT(std::find(tri_conn.begin(), tri_conn.end(), v) != tri_conn.end(),
                             "The vertices of the edge do not belong to the triangle. "
                             "(tri_to_edges_map is wrong)");
            }

            // The position of the edge in the tri_to_edges_map matters, so check that
            // the position follows the order of edges assumed in mesh::get_edge_info.
            _UTIL_ASSERT(
              hash_linear_edge(tri_conn[edge_info[i].connectivity[0]],
                               tri_conn[edge_info[i].connectivity[1]]) ==
                hash_linear_edge(edges[edge_id][0], edges[edge_id][1]),
              "The edge is in the wrong position of the map. (tri_to_edges_map is wrong)");

            // Ensure the triangle appears in edge_to_tris_map for this edge
            bool tri_found = false;
            for (auto t : edge_to_tris_map[edge_id])
            {
                if (t.first == tri_id)
                {
                    tri_found = true;
                }
            }
            _UTIL_ASSERT(tri_found, "The triangle does not appear in edge_to_tris_map.");
        }
    }

    return true;
}

bool check_tri_edge_relations(const Triangulation& t)
{
    return check_tri_edge_relations_kernel(t);
}

bool check_tri_edge_relations(const Tetrahedralization& t)
{
    return check_tri_edge_relations_kernel(t);
}

bool check_tet_edge_relations(const Tetrahedralization& t)
{
    auto& edges            = t.edges;
    auto& tetrahedra       = t.tetrahedra;
    auto& tet_to_edges_map = t.tet_to_edges_map;
    auto& edge_to_tets_map = t.edge_to_tets_map;

    // Verify tetrahedra-edge relationship maps
    for (auto tet_id : utility::range(tetrahedra.size()))
    {
        auto& edge_info = mesh::get_edge_info(mesh::ElementType::T3D4N);
        assert(edge_info.size() == 6);

        // Loop over edges in the tetrahedron
        for (auto i : utility::range(edge_info.size()))
        {
            assert(edge_info[i].type == mesh::EdgeType::E2N);

            auto tet_conn = tetrahedra[tet_id];
            auto& edge_id = tet_to_edges_map[tet_id][i];

            // Check that every vertex in the edge exists in the triangle
            for (auto v : edges[edge_id])
            {
                _UTIL_ASSERT(std::find(tet_conn.begin(), tet_conn.end(), v) != tet_conn.end(),
                             "The vertices of the edge do not belong to the tetrahedron. "
                             "(tet_to_edges_map is wrong)");
            }

            // The position of the edge in the tri_to_edges_map matters, so check that
            // the position follows the order of edges assumed in mesh::get_edge_info.
            _UTIL_ASSERT(
              hash_linear_edge(tet_conn[edge_info[i].connectivity[0]],
                               tet_conn[edge_info[i].connectivity[1]]) ==
                hash_linear_edge(edges[edge_id][0], edges[edge_id][1]),
              "The edge is in the wrong position of the map. (tet_to_edges_map is wrong)");

            // Ensure the triangle appears in edge_to_tris_map for this edge
            bool tet_found = false;
            for (auto t : edge_to_tets_map[edge_id])
            {
                if (t.first == tet_id)
                {
                    tet_found = true;
                }
            }
            if (!tet_found)
                _UTIL_ASSERT(tet_found, "The tetrahedron does not appear in edge_to_tets_map.");
        }
    }

    return true;
}

bool check_tet_tri_relations(const Tetrahedralization& t)
{
    auto& triangles      = t.triangles;
    auto& tetrahedra     = t.tetrahedra;
    auto& tet_to_tri_map = t.tet_to_tri_map;
    auto& tri_to_tet_map = t.tri_to_tet_map;

    // Verify tetrahedra-edge relationship maps
    for (auto tet_id : utility::range(tetrahedra.size()))
    {
        auto& face_info = mesh::get_face_info(mesh::ElementType::T3D4N);
        assert(face_info.size() == 4);

        // Loop over triangles in the tetrahedron
        for (auto i : utility::range(face_info.size()))
        {
            assert(face_info[i].type == mesh::FaceType::T3D3N);

            auto& tet_conn = tetrahedra[tet_id];
            auto& tri_id   = tet_to_tri_map[tet_id][i];

            if (tri_id == mesh::NOT_SET_LocalIndex)
            {
                continue;
            }

            // Check that every vertex in the edge exists in the triangle
            for (auto v : triangles[tri_id])
            {
                _UTIL_ASSERT(std::find(tet_conn.begin(), tet_conn.end(), v) != tet_conn.end(),
                             "The vertices of the triangle do not belong to the tetrahedron. "
                             "(tet_to_tri_map is wrong)");
            }

            // Check that sorted connectivity lists match
            std::array<mesh::LocalIndex, 3> sorted_tri_conn_1 = {
                tet_conn[face_info[i].connectivity[0]],
                tet_conn[face_info[i].connectivity[1]],
                tet_conn[face_info[i].connectivity[2]]
            };
            std::sort(sorted_tri_conn_1.begin(), sorted_tri_conn_1.end());
            auto tri_conn                                     = triangles[tri_id];
            std::array<mesh::LocalIndex, 3> sorted_tri_conn_2 = { tri_conn[0],
                                                                  tri_conn[1],
                                                                  tri_conn[2] };
            std::sort(sorted_tri_conn_2.begin(), sorted_tri_conn_2.end());

            // The position of the edge in the tri_to_edges_map matters, so check that
            // the position follows the order of edges assumed in mesh::get_edge_info.
            _UTIL_ASSERT(
              sorted_tri_conn_1 == sorted_tri_conn_2,
              "The triangle is in the wrong position of the map. (tet_to_tri_map is wrong)");

            // Ensure the tetrahedron appears in tri_to_tet_map for this triangle
            _UTIL_ASSERT(tri_to_tet_map[tri_id][0] == tet_id || tri_to_tet_map[tri_id][1] == tet_id,
                         "The triangle does not appear in tri_to_tet_map.");
        }
    }

    for (auto tri_id : utility::range(triangles.size()))
    {
        auto& tri_tets = tri_to_tet_map[tri_id];
        for (auto& tet_id : tri_tets)
        {
            if (tet_id == mesh::NOT_SET_LocalIndex)
            {
                continue;
            }

            auto& tri_conn = triangles[tri_id];
            auto& tet_conn = tetrahedra[tet_id];

            // Check that every vertex in the edge exists in the triangle
            for (auto v : tri_conn)
            {
                _UTIL_ASSERT(std::find(tet_conn.begin(), tet_conn.end(), v) != tet_conn.end(),
                             "The vertices of the triangle do not belong to the tetrahedron. "
                             "(tri_to_tet_map is wrong)");
            }

            // Ensure the tetrahedron appears in tri_to_tet_map for this triangle
            _UTIL_ASSERT(std::find(tet_to_tri_map[tet_id].begin(),
                                   tet_to_tri_map[tet_id].end(),
                                   tri_id) != tet_to_tri_map[tet_id].end(),
                         "The triangle does not appear in tet_to_tri_map.");
        }
    }

    return true;
}
}