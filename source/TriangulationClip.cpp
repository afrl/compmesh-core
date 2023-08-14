#include "comp_mesh_core/TriangulationClip.h"

#include "comp_mesh_core/SplitEdges.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Polygons.h"
#include "simple_math/Utility.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <deque>
#include <limits>

namespace composite_mesher
{
template<typename V>
Eigen::Vector<boost::multiprecision::cpp_dec_float_50, 3> to_boost_vec(const V& v)
{
    assert(v.size() == 3);
    return Eigen::Vector<boost::multiprecision::cpp_dec_float_50, 3>(v(0), v(1), v(2));
}

namespace plane_line_geom_detail
{
enum class IntersectionType
{
    AT_S1,
    INTERIOR,
    AT_S2,
    LINE_PARALLEL_ON_PLANE,
    LINE_PARALLEL_OFF_PLANE
};

struct IntersectionResult
{
    IntersectionType type;
    simple_math::Point point = { 0., 0., 0. };
};

template<typename V>
IntersectionResult intersect_plane_line_segment_mp(const V& s1_d,
                                                   const V& s2_d,
                                                   const V& p_d,
                                                   const V& n_d)
{
    using boost_float = boost::multiprecision::cpp_dec_float_50;

    auto s1 = to_boost_vec(s1_d);
    auto s2 = to_boost_vec(s2_d);
    auto p  = to_boost_vec(p_d);
    auto n  = to_boost_vec(n_d);

    auto n_dot_s = n.dot(s2 - s1);
    if (n_dot_s.is_zero())
    { // s is parallel to n
        if (n.dot(s1) == n.dot(p))
        {
            return { IntersectionType::LINE_PARALLEL_ON_PLANE };
        }
        else
        {
            return { IntersectionType::LINE_PARALLEL_OFF_PLANE };
        }
    }

    boost_float t = n.dot(p - s1) / n_dot_s;
    if (t.is_zero())
    {
        return { IntersectionType::AT_S1 };
    }
    else if (boost_float(t - boost_float(1)).is_zero())
    {
        return { IntersectionType::AT_S2 };
    }
    else
    {
        Eigen::Vector<boost_float, 3> intersect = s1 + t * (s2 - s1);
        return { IntersectionType::INTERIOR,
                 { intersect(0).convert_to<double>(),
                   intersect(1).convert_to<double>(),
                   intersect(2).convert_to<double>() } };
    }
}

template<typename V>
IntersectionResult intersect_plane_line_segment(const V& s1, const V& s2, const V& p, const V& n)
{
    // Helpful figure to visualize the operations:
    //
    // o s2
    //  \ edge s
    //   \ 
    //    o s1     ^ n, plane_normal
    //             |
    // ------------o----------  plane
    //             p, plane_point

    constexpr auto machine_eps = std::numeric_limits<double>::epsilon();

    auto n_dot_s = n.dot(s2 - s1);
    if (fabs(n_dot_s) < machine_eps)
    {
        // s might be parallel to n, but round off error might be an issue
        return intersect_plane_line_segment_mp(s1, s2, p, n);
    }

    auto t = n.dot(p - s1) / n_dot_s;
    if (t == 0.0)
    {
        return { IntersectionType::AT_S1 };
    }
    else if (fabs(t) < machine_eps)
    {
        // Round off error might be an issue
        return intersect_plane_line_segment_mp(s1, s2, p, n);
    }
    else if (t == 1.0)
    {
        return { IntersectionType::AT_S2 };
    }
    else if (fabs(t - 1.0) < machine_eps)
    {
        // Round off error might be an issue
        return intersect_plane_line_segment_mp(s1, s2, p, n);
    }
    else
    {
        simple_math::Vector3D intersect = s1 + t * (s2 - s1);
        return { IntersectionType::INTERIOR, simple_math::to_point(intersect) };
    }
}

}

namespace plane_point_geom_detail
{
template<typename V>
double point_to_plane_signed_dist(const V& point, const V& plane_point, const V& plane_normal)
{
    constexpr auto machine_eps = std::numeric_limits<double>::epsilon();
    simple_math::Vector3D w    = point - plane_point;
    if (w.norm() < machine_eps)
    {
        using boost_float = boost::multiprecision::cpp_dec_float_50;
        using Vector3DMP  = Eigen::Vector<boost_float, 3>;

        Vector3DMP point_mp(point(0), point(1), point(2));
        Vector3DMP plane_point_mp(plane_point(0), plane_point(1), plane_point(2));
        Vector3DMP plane_normal_mp(plane_normal(0), plane_normal(1), plane_normal(2));

        return plane_normal_mp.dot(point_mp - plane_point_mp).convert_to<double>();
    }
    else
    {
        return plane_normal.dot(w);
    }
}
}

namespace triangle_clip_tables
{
struct Vertex
{
    enum VertexType
    {
        ORIGINAL_TRI_VERTEX,
        NEW_SPLIT_EDGE_VERTEX
    };

    // Either the original vertex index or edge index that is split (depending on type)
    Index i         = NOT_SET_Index;
    VertexType type = ORIGINAL_TRI_VERTEX;
};

struct LookupEntry
{
    // Connectivity of first positive triangle
    std::array<Vertex, 3> pos_tri_0;
    // Connectivity of second positive triangle (often not present)
    std::array<Vertex, 3> pos_tri_1;
    // Connectivity of first negative triangle
    std::array<Vertex, 3> neg_tri_0;
    // Connectivity of second negative triangle (often not present)
    std::array<Vertex, 3> neg_tri_1;
    // Edge aligned with triangle normal x plane normal
    std::array<Vertex, 2> aligned_edge;
};

// Note: 0b111001 = (int) 57, but not all cases are possible.
using LookupTable = std::array<LookupEntry, 0b111001>;

inline auto get_lookup_table_kernel()
{
    LookupTable table;

    constexpr Vertex v0 = { 0, Vertex::ORIGINAL_TRI_VERTEX };
    constexpr Vertex v1 = { 1, Vertex::ORIGINAL_TRI_VERTEX };
    constexpr Vertex v2 = { 2, Vertex::ORIGINAL_TRI_VERTEX };
    constexpr Vertex e0 = { 0, Vertex::NEW_SPLIT_EDGE_VERTEX };
    constexpr Vertex e1 = { 1, Vertex::NEW_SPLIT_EDGE_VERTEX };
    constexpr Vertex e2 = { 2, Vertex::NEW_SPLIT_EDGE_VERTEX };

    table[0b001] = { { v0, e0, e2 }, {}, { e0, v2, e2 }, { e0, v1, v2 }, { e0, e2 } };
    table[0b010] = { { e0, v1, e1 }, {}, { v0, e0, e1 }, { v0, e1, v2 }, { e1, e0 } };
    table[0b100] = { { e1, v2, e2 }, {}, { v0, v1, e1 }, { v0, e1, e2 }, { e2, e1 } };
    table[0b011] = { { v0, v1, e1 }, { v0, e1, e2 }, { e1, v2, e2 }, {}, { e1, e2 } };
    table[0b101] = { { v0, e0, e1 }, { v0, e1, v2 }, { e0, v1, e1 }, {}, { e0, e1 } };
    table[0b110] = { { e0, v2, e2 }, { e0, v1, v2 }, { v0, e0, e2 }, {}, { e2, e0 } };

    table[0b010001] = { { v0, v1, e2 }, {}, { v1, v2, e2 }, {}, { v1, e2 } };
    table[0b100001] = { { v0, e0, v2 }, {}, { e0, v1, v2 }, {}, { e0, v2 } };
    table[0b001010] = { { v0, v1, e1 }, {}, { v0, e1, v2 }, {}, { e1, v0 } };
    table[0b100010] = { { e0, v1, v2 }, {}, { v0, e0, v2 }, {}, { v2, e0 } };
    table[0b001100] = { { v0, e1, v2 }, {}, { v0, v1, e1 }, {}, { v0, e1 } };
    table[0b010100] = { { e2, v1, v2 }, {}, { v0, v1, e2 }, {}, { e2, v1 } };

    table[0b110001] = { { v0, v1, v2 }, {}, {}, {}, { e1, e2 } };
    table[0b110000] = { {}, { v0, v1, v2 }, {}, {}, { e2, e1 } };
    table[0b101010] = { { v0, v1, v2 }, {}, {}, {}, { e2, e0 } };
    table[0b101000] = { {}, { v0, v1, v2 }, {}, {}, { e0, e2 } };
    table[0b011100] = { { v0, v1, v2 }, {}, {}, {}, { e0, e1 } };
    table[0b011000] = { {}, { v0, v1, v2 }, {}, {}, { e1, e0 } };

    return table;
}

const auto& get_lookup_table()
{
    static auto table = get_lookup_table_kernel();
    return table;
}
}

// New algorithm:
// 1. Compute signed distance field at vertices (use multiprecision)
// 2. Loop over triangles
//    - Skip those that are all negative
//    - Add those that are all positive
//    - Use table to split thse that have changes in sign (track new and existing edges on plane)
//      Note: reorient edges on plane to point outward using triangle information since triangles
//            are already oriented.
//    - Add subtriangles that are all positive
// 3. Find polyline loops (ordered and oriented)
// 4. Triangulate on plane and add

///////////////////////////////////////////////////////////////////////////////////////////////
// Definitions from header

using namespace plane_line_geom_detail;
using namespace plane_point_geom_detail;

Triangulation clip_by_plane(const Triangulation& t,
                            const simple_math::Point& plane_point,
                            const simple_math::Point& plane_normal,
                            const bool& close_clipped_surfaces)
{
    Triangulation clipped_t;

    // Compute signed distance field at vertices
    auto p = simple_math::point_as_vector(plane_point);
    auto n = simple_math::point_as_vector(plane_normal);
    std::vector<double> vert_signed_dist(t.vertices.size(), 0.0);
    std::vector<Index> old_to_new_v_map(t.vertices.size(), NOT_SET_Index);
    for (auto v_id : utility::indices(t.vertices))
    {
        auto v                 = simple_math::point_as_vector(t.vertices[v_id]);
        vert_signed_dist[v_id] = point_to_plane_signed_dist(v, p, n);
        if (vert_signed_dist[v_id] <= 0.0)
        {
            old_to_new_v_map[v_id] = clipped_t.vertices.size();
            clipped_t.vertices.push_back(t.vertices[v_id]);
        }
    }

    // For the edges that are intersected by the plane, compute the intersections
    std::map<Index, Index> edge_id_to_intersect_v_id_map;
    for (auto edge_id : utility::indices(t.edges))
    {
        const auto& edge = t.edges[edge_id];
        if (simple_math::sign(vert_signed_dist[edge[0]]) !=
            simple_math::sign(vert_signed_dist[edge[1]]))
        {
            auto s1 = simple_math::point_as_vector(t.vertices[edge[0]]);
            auto s2 = simple_math::point_as_vector(t.vertices[edge[1]]);
            auto r  = intersect_plane_line_segment(s1, s2, p, n);
            assert(r.type == IntersectionType::INTERIOR);
            edge_id_to_intersect_v_id_map[edge_id] = clipped_t.vertices.size();
            clipped_t.vertices.push_back(r.point);
        }
    }

    const auto& lookup_table = triangle_clip_tables::get_lookup_table();
    std::vector<std::array<Index, 2>> edges_on_plane;
    for (auto tri_id : utility::indices(t.triangles))
    {
        auto& tri = t.triangles[tri_id];
        if (std::all_of(tri.begin(),
                        tri.end(),
                        [&vert_signed_dist](auto& v_id) -> bool
                        { return vert_signed_dist[v_id] < 0.0; }))
        {
            // Triangle is on the negative side of the plane, so keep it.
            std::array<Index, 3> new_tri;
            for (auto j : utility::range(3))
            {
                assert(old_to_new_v_map[tri[j]] != NOT_SET_Index);
                new_tri[j] = old_to_new_v_map[tri[j]];
            }
            clipped_t.triangles.push_back(new_tri);
            clipped_t.tri_to_marker_map.push_back(t.tri_to_marker_map[tri_id]);
        }
        else if (std::all_of(tri.begin(),
                             tri.end(),
                             [&vert_signed_dist](auto& v_id) -> bool
                             { return vert_signed_dist[v_id] == 0.0; }))
        {
            // Triangle is on the plane...
            throw _UTIL_MSG_EXCEPT("TODO implement");
        }
        else if (!std::all_of(tri.begin(),
                              tri.end(),
                              [&vert_signed_dist](auto& v_id) -> bool
                              { return vert_signed_dist[v_id] > 0.0; }))
        {
            // Triangle is cut by the plane
            // The table key is effectively made of two bitsets:
            // |__|__|__|__|__|__|
            //   node i | node i
            //   == iso | > iso
            uint8_t key = 0;
            // > iso section
            for (auto v : utility::range(3))
            {
                const auto& v_id   = tri[v];
                const uint8_t flag = vert_signed_dist[v_id] > 0.0;
                key |= flag << v;
            }
            // == iso section
            for (auto v : utility::range(3))
            {
                const auto& v_id   = tri[v];
                const uint8_t flag = vert_signed_dist[v_id] == 0.0;
                key |= flag << (v + 3);
            }
            const auto& entry = lookup_table[key];

            std::array<Index, 3> new_tri;
            for (const auto& neg_tri : { entry.neg_tri_0, entry.neg_tri_1 })
            {
                // Keep the triangles on the negative side of the plane
                if (neg_tri[0].i != NOT_SET_Index)
                {
                    for (auto j : utility::range(3))
                    {
                        const auto& index = neg_tri[j].i;
                        const auto& type  = neg_tri[j].type;

                        if (type == neg_tri[j].ORIGINAL_TRI_VERTEX)
                        {
                            assert(old_to_new_v_map[tri[index]] != NOT_SET_Index);
                            new_tri[j] = old_to_new_v_map[tri[index]];
                        }
                        else
                        {
                            auto find_iter =
                              edge_id_to_intersect_v_id_map.find(t.tri_to_edges_map[tri_id][index]);
                            assert(find_iter != edge_id_to_intersect_v_id_map.end());
                            new_tri[j] = find_iter->second;
                        }
                    }
                    clipped_t.triangles.push_back(new_tri);
                    clipped_t.tri_to_marker_map.push_back(t.tri_to_marker_map[tri_id]);
                }
            }

            std::array<Index, 2> new_edge;
            for (auto j : utility::range(2))
            {
                const auto& index = entry.aligned_edge[j].i;
                const auto& type  = entry.aligned_edge[j].type;

                if (type == entry.aligned_edge[j].ORIGINAL_TRI_VERTEX)
                {
                    assert(old_to_new_v_map[tri[index]] != NOT_SET_Index);
                    new_edge[j] = old_to_new_v_map[tri[index]];
                }
                else
                {
                    auto find_iter =
                      edge_id_to_intersect_v_id_map.find(t.tri_to_edges_map[tri_id][index]);
                    assert(find_iter != edge_id_to_intersect_v_id_map.end());
                    new_edge[j] = find_iter->second;
                }
            }
            edges_on_plane.push_back(new_edge);
        }
        // else: triangle is on the positive side of the plane, so skip it.
    }

    // Find the polygons created by the clip operation and triangulate them
    if (close_clipped_surfaces)
    {
        auto insert_into_map =
          [](std::map<Index, std::array<Index, 2>>& m, const Index& v_id, const Index& edge_id)
        {
            auto find_iter = m.find(v_id);
            if (find_iter == m.end())
            {
                m[v_id] = { edge_id, NOT_SET_Index };
            }
            else
            {
                assert(find_iter->second[1] == NOT_SET_Index);
                find_iter->second[1] = edge_id;
            }
        };

        std::map<Index, std::array<Index, 2>> vertex_id_to_edge_on_plane_map;
        for (Index edge_id : utility::indices(edges_on_plane))
        {
            auto& edge     = edges_on_plane[edge_id];
            auto find_iter = vertex_id_to_edge_on_plane_map.find(edge[0]);
            if (find_iter == vertex_id_to_edge_on_plane_map.end())
            {
                vertex_id_to_edge_on_plane_map[edge[0]] = { edge_id, NOT_SET_Index };
            }
            else
            {
                assert(find_iter->second[0] == NOT_SET_Index);
                find_iter->second[0] = edge_id;
            }
            find_iter = vertex_id_to_edge_on_plane_map.find(edge[1]);
            if (find_iter == vertex_id_to_edge_on_plane_map.end())
            {
                vertex_id_to_edge_on_plane_map[edge[1]] = { NOT_SET_Index, edge_id };
            }
            else
            {
                assert(find_iter->second[1] == NOT_SET_Index);
                find_iter->second[1] = edge_id;
            }
        }

        // Traverse edges_on_plane to find loops
        Index n_edges_added     = 0;
        Index next_start_edge_i = 0;
        std::vector<bool> is_edge_on_plane_added(edges_on_plane.size(), false);
        std::vector<std::vector<Index>> vertex_loops;
        while (next_start_edge_i != edges_on_plane.size())
        {
            std::vector<Index> curr_loop              = { edges_on_plane[next_start_edge_i][0],
                                                          edges_on_plane[next_start_edge_i][1] };
            is_edge_on_plane_added[next_start_edge_i] = true;

            while (true)
            {
                auto& next_edge_i = vertex_id_to_edge_on_plane_map[curr_loop.back()][0];
                is_edge_on_plane_added[next_edge_i] = true;
                auto& next_v_id                     = edges_on_plane[next_edge_i][1];
                if (next_v_id == curr_loop.front())
                {
                    break;
                }
                curr_loop.push_back(next_v_id);
            }
            vertex_loops.push_back(curr_loop);

            for (; next_start_edge_i < edges_on_plane.size(); next_start_edge_i++)
            {
                if (!is_edge_on_plane_added[next_start_edge_i])
                {
                    break;
                }
            }
        }

        // Submesh polygon loops
        {
            /*
            {
                mesh::Mesh m;
                for (auto& v : clipped_t.vertices)
                {
                    m.nodes().add(v);
                }
                for (auto i : utility::indices(vertex_loops))
                {
                    for (auto j : utility::range((std::size_t)1, vertex_loops[i].size()))
                    {
                        m.elements().add(mesh::ElementType::L1D2N,
                                         { vertex_loops[i][j - 1], vertex_loops[i][j] });
                        m.elements().back().component_id() = i;
                    }
                    m.elements().add(mesh::ElementType::L1D2N,
                                     { vertex_loops[i].back(), vertex_loops[i][0] });
                    m.elements().back().component_id() = i;
                }
                mesh::write_vtk(m, "loops.vtk");
            }
            */

            mesh::Mesh m;
            for (auto& v : clipped_t.vertices)
            {
                m.nodes().add(v);
            }
            for (auto i : utility::indices(vertex_loops))
            {

                m.elements().add(mesh::ElementType::POLYGON, vertex_loops[i]);
            }

            // mesh::write_vtk(m, "caps.vtk");

            triangulate_all_polygons(m, -1);

            // mesh::write_vtk(m, "caps_post_triangulate.vtk");

            // Note: triangulation may add some points
            for (auto i : utility::range(clipped_t.vertices.size(), (std::size_t)m.nodes().size()))
            {
                clipped_t.vertices.push_back(m.nodes()[i].point());
            }

            for (auto& e : m.elements())
            {
                assert(e.type() == mesh::ElementType::T2D3N);
                auto conn = e.connectivity();
                clipped_t.triangles.push_back({ conn[0], conn[1], conn[2] });
                clipped_t.tri_to_marker_map.push_back(clipped_t.n_markers);
            }
            clipped_t.n_markers++;
        }
    }

    update_edges_and_maps(clipped_t);

    check_tri_edge_relations(clipped_t);

    return clipped_t;
}
}