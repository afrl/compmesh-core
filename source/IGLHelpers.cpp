#include "comp_mesh_core/IGLHelpers.h"

#include "mesh/algorithms/Search.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"
#include "utility/Range.h"

#include <unordered_set>

namespace composite_mesher
{
namespace detail
{
template<typename V, typename M>
IGLToMeshResult igl_to_mesh_kernel(const V& vertices, const M& triangles)
{
    IGLToMeshResult r;
    for (auto i : utility::range(vertices.rows()))
    {
        r.mesh.nodes().add({ vertices(i, 0), vertices(i, 1), vertices(i, 2) });
    }
    r.element_to_origin_tris_map.reserve(triangles.rows());
    std::unordered_map<std::size_t, mesh::LocalIndex> inserted_triangles;
    for (auto i : utility::range(triangles.rows()))
    {
        std::vector<mesh::LocalIndex> conn = { static_cast<mesh::LocalIndex>(triangles(i, 0)),
                                               static_cast<mesh::LocalIndex>(triangles(i, 1)),
                                               static_cast<mesh::LocalIndex>(triangles(i, 2)) };
        auto hash                          = utility::hash_sorted(conn);
        auto find_iter                     = inserted_triangles.find(hash);
        if (find_iter == inserted_triangles.end())
        {
            inserted_triangles[hash] = r.mesh.elements().size();
            r.mesh.elements().add(mesh::ElementType::T2D3N, conn);
            r.element_to_origin_tris_map.push_back({ i });
        }
        else
        {
            r.element_to_origin_tris_map[find_iter->second].push_back(i);
        }
    }
    return r;
}

}
IGLTriSurf mesh_to_igl(const mesh::Mesh& m)
{
    IGLTriSurf surf;
    surf.vertices.resize(m.nodes().size(), 3);
    for (auto& n : m.nodes())
    {
        for (auto i : utility::range(3))
        {
            surf.vertices(n.local_index(), i) = n.point(i);
        }
    }
    surf.triangles.resize(m.elements().size(), 3);
    for (auto& e : m.elements())
    {
        auto& conn = e.connectivity();
        for (auto i : utility::range(3))
        {
            surf.triangles(e.local_index(), i) = conn[i];
        }
    }

    surf.tree.init(surf.vertices, surf.triangles);

    return surf;
}

IGLTriSurf mesh_to_igl(const mesh::Mesh& m, const mesh::ComponentIndex& component_to_extract)
{
    mesh::LocalIndex curr_index              = 0;
    mesh::LocalIndex n_elements_in_component = 0;
    std::vector<mesh::LocalIndex> m_to_c_node_map(m.nodes().size(),
                                                  std::numeric_limits<mesh::LocalIndex>::max());
    std::vector<mesh::LocalIndex> c_to_m_node_map;
    c_to_m_node_map.reserve(m.nodes().size());
    for (auto& e : m.elements())
    {
        if (e.component_id() == component_to_extract)
        {
            n_elements_in_component++;
            auto& conn = e.connectivity();
            assert(conn.size() == 3); // only works for triangles
            for (auto i : utility::range(3))
            {
                if (m_to_c_node_map[conn[i]] == std::numeric_limits<mesh::LocalIndex>::max())
                {
                    m_to_c_node_map[conn[i]] = curr_index++;
                    c_to_m_node_map.push_back(conn[i]);
                }
            }
        }
    }

    IGLTriSurf surf;
    surf.vertices.resize(curr_index, 3);
    for (auto i : utility::indices(c_to_m_node_map))
    {
        auto n = m.nodes()[c_to_m_node_map[i]];
        for (auto j : utility::range(3))
        {
            surf.vertices(i, j) = n.point(j);
        }
    }

    curr_index = 0;
    surf.triangles.resize(n_elements_in_component, 3);
    for (auto& e : m.elements())
    {
        if (e.component_id() == component_to_extract)
        {
            auto& conn = e.connectivity();
            for (auto i : utility::range(3))
            {
                surf.triangles(curr_index, i) = m_to_c_node_map[conn[i]];
            }
            curr_index++;
        }
    }

    surf.tree.init(surf.vertices, surf.triangles);

    return surf;
}

IGLTriSurf mesh_to_igl(const mesh::Mesh& m,
                       const std::vector<mesh::ComponentIndex>& components_to_extract)
{
    std::unordered_set<mesh::ComponentIndex> components_to_extract_set(
      components_to_extract.begin(),
      components_to_extract.end());

    mesh::LocalIndex curr_index              = 0;
    mesh::LocalIndex n_elements_in_component = 0;
    std::vector<mesh::LocalIndex> m_to_c_node_map(m.nodes().size(),
                                                  std::numeric_limits<mesh::LocalIndex>::max());
    std::vector<mesh::LocalIndex> c_to_m_node_map;
    c_to_m_node_map.reserve(m.nodes().size());
    for (auto& e : m.elements())
    {
        if (components_to_extract_set.count(e.component_id()) > 0)
        {
            n_elements_in_component++;
            auto& conn = e.connectivity();
            assert(conn.size() == 3); // only works for triangles
            for (auto i : utility::range(3))
            {
                if (m_to_c_node_map[conn[i]] == std::numeric_limits<mesh::LocalIndex>::max())
                {
                    m_to_c_node_map[conn[i]] = curr_index++;
                    c_to_m_node_map.push_back(conn[i]);
                }
            }
        }
    }

    IGLTriSurf surf;
    surf.vertices.resize(curr_index, 3);
    for (auto i : utility::indices(c_to_m_node_map))
    {
        auto n = m.nodes()[c_to_m_node_map[i]];
        for (auto j : utility::range(3))
        {
            surf.vertices(i, j) = n.point(j);
        }
    }

    curr_index = 0;
    surf.triangles.resize(n_elements_in_component, 3);
    for (auto& e : m.elements())
    {
        if (components_to_extract_set.count(e.component_id()) > 0)
        {
            auto& conn = e.connectivity();
            for (auto i : utility::range(3))
            {
                surf.triangles(curr_index, i) = m_to_c_node_map[conn[i]];
            }
            curr_index++;
        }
    }

    surf.tree.init(surf.vertices, surf.triangles);

    return surf;
}

IGLTriSurf tetrahedralization_to_igl(const Tetrahedralization& t,
                                     const std::vector<mesh::LocalIndex>& triangles_to_extract)
{
    mesh::LocalIndex curr_index = 0;
    std::vector<mesh::LocalIndex> m_to_c_node_map(t.vertices.size(),
                                                  std::numeric_limits<mesh::LocalIndex>::max());
    std::vector<mesh::LocalIndex> c_to_m_node_map;
    c_to_m_node_map.reserve(t.vertices.size());
    for (auto& tri_id : triangles_to_extract)
    {
        auto& conn = t.triangles[tri_id];
        assert(conn.size() == 3); // only works for triangles
        for (auto i : utility::range(3))
        {
            if (m_to_c_node_map[conn[i]] == std::numeric_limits<mesh::LocalIndex>::max())
            {
                m_to_c_node_map[conn[i]] = curr_index++;
                c_to_m_node_map.push_back(conn[i]);
            }
        }
    }


    IGLTriSurf surf;
    surf.vertices.resize(curr_index, 3);
    for (auto i : utility::indices(c_to_m_node_map))
    {
        auto& p = t.vertices[c_to_m_node_map[i]];
        for (auto j : utility::range(3))
        {
            surf.vertices(i, j) = p[j];
        }
    }

    curr_index = 0;
    surf.triangles.resize(triangles_to_extract.size(), 3);
    for (auto& tri_id : triangles_to_extract)
    {
        auto& conn = t.triangles[tri_id];
        for (auto i : utility::range(3))
        {
            surf.triangles(curr_index, i) = m_to_c_node_map[conn[i]];
        }
        curr_index++;
    }

    surf.tree.init(surf.vertices, surf.triangles);

    return surf;
}



IGLToMeshResult igl_to_mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& triangles)
{
    return detail::igl_to_mesh_kernel(vertices, triangles);
}

IGLToMeshResult igl_to_mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixX3i& triangles)
{
    return detail::igl_to_mesh_kernel(vertices, triangles);
}

double igl_sq_distance(const IGLTriSurf& igl_surf, const mesh::Point& p)
{
    Eigen::MatrixX3d point(1, 3);
    point.row(0) = simple_math::point_as_vector(p);

    Eigen::VectorXd squared_distance;
    Eigen::VectorXi closest_triangle;
    Eigen::MatrixX3d closest_point;
    igl_surf.tree.squared_distance(igl_surf.vertices,
                                   igl_surf.triangles,
                                   point,
                                   squared_distance,
                                   closest_triangle,
                                   closest_point);
    return squared_distance(0);
}

Eigen::Vector3d igl_closest_point(const IGLTriSurf& igl_surf, const Eigen::Vector3d& p)
{
    Eigen::MatrixX3d point(1, 3);
    point.row(0) = p;

    Eigen::VectorXd squared_distance;
    Eigen::VectorXi closest_triangle;
    Eigen::MatrixX3d closest_point;
    igl_surf.tree.squared_distance(igl_surf.vertices,
                                   igl_surf.triangles,
                                   point,
                                   squared_distance,
                                   closest_triangle,
                                   closest_point);

    return closest_point.row(0);
}

std::vector<IGLTriSurf> decompose_into_igl_facets(mesh::Mesh combined_surface_mesh)
{
    // Convert each facet of the input surfaces into triangle surfaces in IGL's format
    std::vector<IGLTriSurf> facet_id_to_igl_surface;
    for (auto facet_id : utility::range(combined_surface_mesh.n_components()))
    {
        // Note: I am assuming that the components are repurposed for tracking the facets
        facet_id_to_igl_surface.push_back(mesh_to_igl(combined_surface_mesh, facet_id));
    }
    return facet_id_to_igl_surface;
}

void igl_closest_point_to_facet(const Eigen::MatrixX3d& vertices,
                                const IGLTriSurf& igl_surf,
                                Eigen::MatrixX3d& vectors_to_surface,
                                Eigen::VectorXd& distances_to_surface)
{
    Eigen::VectorXd squared_distances;
    Eigen::VectorXi closest_triangles;
    Eigen::MatrixX3d closest_points;
    igl_surf.tree.squared_distance(igl_surf.vertices,
                                   igl_surf.triangles,
                                   vertices,
                                   squared_distances,
                                   closest_triangles,
                                   closest_points);

    vectors_to_surface = closest_points - vertices;
    // squared_distances.cwiseSqrt();
    distances_to_surface = squared_distances.cwiseSqrt();
}

IGLTetrahedralization mesh_to_igl_tet(const mesh::Mesh& m)
{
    IGLTetrahedralization igl;
    igl.vertices.resize(m.nodes().size(), 3);
    for (auto& n : m.nodes())
    {
        for (auto i : utility::range(3))
        {
            igl.vertices(n.local_index(), i) = n.point(i);
        }
    }
    igl.tetrahedra.resize(m.elements().size(), 4);
    for (auto& e : m.elements())
    {
        auto& conn = e.connectivity();
        for (auto i : utility::range(4))
        {
            igl.tetrahedra(e.local_index(), i) = conn[i];
        }
    }
    return igl;
}
}