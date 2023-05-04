#include "comp_mesh_core/TriangulationClip.h"

#include "comp_mesh_core/SplitEdges.h"
#include "simple_math/Utility.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <limits>

namespace composite_mesher
{
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

    Eigen::Vector<boost_float, 3> s1 = s1_d;
    Eigen::Vector<boost_float, 3> s2 = s2_d;
    Eigen::Vector<boost_float, 3> p  = p_d;
    Eigen::Vector<boost_float, 3> n  = n_d;

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

        Vector3DMP point_mp        = point;
        Vector3DMP plane_point_mp  = plane_point;
        Vector3DMP plane_normal_mp = plane_normal;

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
struct LookupEntry
{
    enum VertexType
    {
        ORIGINAL_TRI_VERTEX,
        NEW_SPLIT_EDGE_VERTEX
    };

    // Connectivity of first positive triangle
    std::array<Index, 3> pos_tri_0 = { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index };
    std::array<VertexType, 3> pos_tri_0_vert_types;
    // Connectivity of second positive triangle (often not present)
    std::array<Index, 3> pos_tri_1 = { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index };
    std::array<VertexType, 3> pos_tri_1_vert_types;
    // Connectivity of first negative triangle
    std::array<Index, 3> neg_tri_0 = { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index };
    std::array<VertexType, 3> neg_tri_0_vert_types;
    // Connectivity of second negative triangle (often not present)
    std::array<Index, 3> neg_tri_1 = { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index };
    std::array<VertexType, 3> neg_tri_1_vert_types;
    // Edge aligned with triangle normal x plane normal
    std::array<Index, 2> aligned_edge = { NOT_SET_Index, NOT_SET_Index };
};

// Note: 0b111001 = (int) 57, but not all cases are possible.
using LookupTable = std::array<LookupEntry, 0b111001>;

inline auto get_cell_lookup_table_kernel()
{
    LookupTable table;

    table[0b001] = { { 0, e0, e2 } };
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

void clip_by_plane(const Triangulation& t,
                   const simple_math::Point& plane_point,
                   const simple_math::Point& plane_normal)
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
        else if (!std::all_of(tri.begin(),
                              tri.end(),
                              [&vert_signed_dist](auto& v_id) -> bool
                              { return vert_signed_dist[v_id] > 0.0; }))
        {
            // Triangle is cut by the plane
            // TODO
        }
        // else: triangle is on the positive side of the plane, so skip it.
    }



    constexpr auto machine_eps = std::numeric_limits<double>::epsilon();
    auto p                     = simple_math::point_as_vector(plane_point);
    auto n                     = simple_math::point_as_vector(plane_normal);

    std::vector<Index> edges_on_plane;
    std::vector<Index> edges_to_split;
    std::vector<Index> intersection_point_vertices;
    for (auto edge_id : utility::indices(t.edges))
    {
        const auto& edge = t.edges[edge_id];
        auto s1          = simple_math::point_as_vector(t.vertices[edge[0]]);
        auto s2          = simple_math::point_as_vector(t.vertices[edge[1]]);

        auto r = intersect_plane_line_segment(s1, s2, p, n);
        if (r.type == IntersectionType::INTERIOR)
        {
            edges_to_split.push_back(edge_id);
            intersection_point_vertices.push_back(t.vertices.size());
            t.vertices.push_back(r.point);
        }
        else if (r.type == IntersectionType::LINE_PARALLEL_ON_PLANE)
        {
            edges_on_plane.push_back(edge_id);
        }
    }

    Index n_edges_found_on_plane = 0;
    auto update                  = split_edges(edges_to_split, intersection_point_vertices, t);
    for (auto& edge_update : update.updated_edges)
    {
        if (edge_update.update_type == edge_update.NEW_EDGE)
        {

            const auto& edge = t.edges[edge_update.new_edge_id];
            auto s1          = simple_math::point_as_vector(t.vertices[edge[0]]);
            auto s2          = simple_math::point_as_vector(t.vertices[edge[1]]);
            auto r           = intersect_plane_line_segment(s1, s2, p, n);
            if (r.type == IntersectionType::LINE_PARALLEL_ON_PLANE)
            {
                edges_on_plane.push_back(edge_update.new_edge_id);
                n_edges_found_on_plane++;
            }
        }
    }
    assert(n_edges_found_on_plane == edges_to_split.size());
}
}