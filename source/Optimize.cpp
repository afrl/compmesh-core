#include "comp_mesh_core/Optimize.h"

#include "comp_mesh_core/GlobalOptions.h"
#include "comp_mesh_core/IGLHelpers.h"
#include "comp_mesh_core/VertexMap.h"
#include "mesh/algorithms/IO.h"
#include "mesh/relations/NodeToElementMap.h"
#include "simple_math/LineGeometry.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"
#include "utility/Permutation.h"
#include "utility/VectorSet.h"

namespace composite_mesher
{
namespace detail
{

void write_seams(const Triangulation& t,
                 const SeamInfo& seam_info,
                 const std::string prefix  = "",
                 const std::string postfix = "")
{
    mesh::Mesh seams;
    for (auto& p : t.vertices)
    {
        seams.nodes().add(p);
    }

    std::vector<Marker> seam_id;
    for (auto edge_id : seam_info.edges_on_seams)
    {
        seams.elements().add(mesh::ElementType::L1D2N,
                             { t.edges[edge_id][0], t.edges[edge_id][1] });
        seam_id.push_back(seam_info.edge_id_to_seam_id_map[edge_id]);
    }

    mesh::MeshAssociatedData data;
    data.add<mesh::DataFieldType::ELEMENT, Marker>(
      mesh::DataField<Marker>::ref("seam_id", seam_id));

    mesh::write_vtk(seams, prefix + "seams" + postfix + ".vtk", data);
}

void write_seams(const mesh::Mesh& surface_mesh,
                 const std::vector<std::array<Index, 2>>& edges,
                 const Marker& n_seams,
                 const std::vector<Marker>& edges_on_seams,
                 const std::vector<Marker>& edge_id_to_seam_id_map,
                 const std::string prefix  = "",
                 const std::string postfix = "")
{
    mesh::Mesh seams;
    for (auto& n : surface_mesh.nodes())
    {
        seams.nodes().add(n.point());
    }

    std::vector<Marker> seam_id;
    for (auto edge_id : edges_on_seams)
    {
        seams.elements().add(mesh::ElementType::L1D2N, { edges[edge_id][0], edges[edge_id][1] });
        seam_id.push_back(edge_id_to_seam_id_map[edge_id]);
    }

    mesh::MeshAssociatedData data;
    data.add<mesh::DataFieldType::ELEMENT, Marker>(
      mesh::DataField<Marker>::ref("seam_id", seam_id));

    mesh::write_vtk(seams, prefix + "seams" + postfix + ".vtk", data);
}

struct AdjacentNode
{
    Index node_id;
    Index edge_id;
};
}



void optimize_tangential_smoothing(mesh::Mesh& surface_mesh,
                                   const std::vector<std::array<mesh::LocalIndex, 2>>& edges,
                                   const mesh::LocalIndex& n_seams,
                                   const std::vector<mesh::LocalIndex>& edges_on_seams,
                                   const std::vector<mesh::LocalIndex>& edge_id_to_seam_id_map,
                                   const bool allow_smoothing_on_seams)
{
    std::cout << GlobalOptions::get().curr_status_padding
              << "* optimize_tangential_smoothing - Starting iterative optimization" << std::endl;

    if (GlobalOptions::get().write_debug_meshes)
    {
        mesh::write_vtk(surface_mesh, GlobalOptions::get().output_dir + "/iteration_0.vtk");
        detail::write_seams(surface_mesh,
                            edges,
                            n_seams,
                            edges_on_seams,
                            edge_id_to_seam_id_map,
                            GlobalOptions::get().output_dir + "/iteration_0_");
    }

    // Count # of facets and seams
    mesh::ComponentIndex n_facets = surface_mesh.n_components();

    mesh::NodeToElementMap node_to_element_map;
    node_to_element_map.update(surface_mesh);

    // Make a copy of the original node positions
    std::vector<simple_math::Point> original_positions(surface_mesh.nodes().size());
    for (auto n : surface_mesh.nodes())
    {
        original_positions[n.local_index()] = n.point();
    }

    // Create a map from each feature to the nodes that compose them
    std::vector<std::vector<mesh::LocalIndex>> facet_id_to_node_id_map(n_facets);
    for (auto n_id : utility::range(surface_mesh.nodes().size()))
    {
        for (auto e : node_to_element_map.get_elements(n_id))
        {
            facet_id_to_node_id_map[e.component_id()].push_back(n_id);
        }
    }
    std::vector<bool> is_node_on_seam(surface_mesh.nodes().size(), false);
    std::vector<mesh::ComponentIndex> n_seam_elements_connected_to_node(surface_mesh.nodes().size(),
                                                                        0);
    std::vector<std::vector<mesh::LocalIndex>> seam_id_to_node_id_map(n_seams);
    std::vector<std::vector<detail::AdjacentNode>> node_to_neighbors_on_seam(
      surface_mesh.nodes().size());
    for (auto feature_edge_id : edges_on_seams)
    {
        auto seam_id = edge_id_to_seam_id_map[feature_edge_id];

        auto n1_id             = edges[feature_edge_id][0];
        is_node_on_seam[n1_id] = true;
        n_seam_elements_connected_to_node[n1_id]++;
        seam_id_to_node_id_map[seam_id].push_back(n1_id);

        auto n2_id             = edges[feature_edge_id][1];
        is_node_on_seam[n2_id] = true;
        n_seam_elements_connected_to_node[n2_id]++;
        seam_id_to_node_id_map[seam_id].push_back(n2_id);

        node_to_neighbors_on_seam[n1_id].push_back({ n2_id, feature_edge_id });
        node_to_neighbors_on_seam[n2_id].push_back({ n1_id, feature_edge_id });
    }

    // Create a map from each node to the facets it should lie on
    std::vector<utility::VectorSet<unsigned int>> node_id_to_input_facet_ids_map(
      surface_mesh.nodes().size());
    for (auto e_id : utility::range(surface_mesh.elements().size()))
    {
        for (auto n_id : surface_mesh.elements()[e_id].connectivity())
        {
            node_id_to_input_facet_ids_map[n_id].insert(
              surface_mesh.elements()[e_id].component_id());
        }
    }

    // Figure out which nodes are on a corner based on how many seam elements are connected to
    // the node.  It is assumed that if it is greater than 2, then it a corner.
    std::vector<bool> is_node_on_corner(surface_mesh.nodes().size(), false);
    for (auto n_id : utility::range(surface_mesh.nodes().size()))
    {
        if (n_seam_elements_connected_to_node[n_id] > 2)
        {
            is_node_on_corner[n_id] = true;
        }
    }

    // TODO build adjacency list for seam edges

    auto igl_facets = decompose_into_igl_facets(surface_mesh);

    std::vector<uint32_t> is_node_on_seam_data(is_node_on_seam.size(), 0);
    for (auto i : utility::indices(is_node_on_seam_data))
    {
        if (is_node_on_seam[i])
        {
            is_node_on_seam_data[i] = is_node_on_seam[i];
        }
    }

    mesh::MeshAssociatedData data;
    data.add<mesh::DataFieldType::NODE, uint32_t>(
      mesh::DataField<uint32_t>::cref("is_on_seam", is_node_on_seam_data));

    for (auto iteration : utility::range(10))
    {
        std::cout << GlobalOptions::get().curr_status_padding << " - Starting iteration "
                  << iteration << std::endl;

        std::vector<simple_math::Vector3D> updated_coords(surface_mesh.nodes().size(),
                                                          simple_math::Vector3D::Zero());

        if (allow_smoothing_on_seams)
        {
            for (auto seam_id : utility::range(n_seams))
            {
                for (auto& n_id : seam_id_to_node_id_map[seam_id])
                {
                    if (!is_node_on_corner[n_id])
                    {
                        simple_math::Vector3D u = simple_math::Vector3D::Zero();
                        simple_math::Vector3D t = simple_math::Vector3D::Zero();
                        double edge_length_sums = 0;
                        int n_connected_edges   = 0;
                        auto p1                 = surface_mesh.nodes()[n_id].coord();

                        // Loop over points on feature connected to n_id via an edge
                        for (auto neighbor : node_to_neighbors_on_seam[n_id])
                        {
                            auto p2 = surface_mesh.nodes()[neighbor.node_id].coord();

                            simple_math::Vector3D midpoint = 0.5 * (p1 + p2);
                            double distance                = (p2 - p1).norm();
                            edge_length_sums += distance;
                            u += distance * midpoint;

                            if (n_connected_edges == 0)
                            {
                                t += (p2 - p1).normalized();
                                n_connected_edges++;
                            }
                            else
                            {
                                n_connected_edges++;
                                t -= (p2 - p1).normalized();
                            }
                        }

                        assert(n_connected_edges == 2);

                        if (edge_length_sums > 0)
                        { // to avoid overflow (i.e., edge_length_sums == 0)
                            u /= edge_length_sums;
                            u -= p1;
                            t                    = t.normalized();
                            u                    = t * u.dot(t);
                            updated_coords[n_id] = u;
                        }
                    }
                }
            }
        }

        for (auto facet_id : utility::range(n_facets))
        {
            for (auto& n_id : facet_id_to_node_id_map[facet_id])
            {
                if (!is_node_on_corner[n_id] && !is_node_on_seam[n_id])
                {
                    simple_math::Vector3D u = simple_math::Vector3D::Zero();
                    simple_math::Vector3D t = simple_math::Vector3D::Zero();
                    double ww               = 0;
                    int c                   = 0;
                    auto p1                 = surface_mesh.nodes()[n_id].coord();

                    for (auto parent_element : node_to_element_map.get_elements(n_id))
                    {
                        if (parent_element.type() == mesh::ElementType::T2D3N &&
                            parent_element.component_id() == facet_id)
                        {
                            auto p1   = parent_element.node(0).coord();
                            auto p2   = parent_element.node(1).coord();
                            auto p3   = parent_element.node(2).coord();
                            auto avg  = (p1 + p2 + p3) / 3.0;
                            auto area = (p2 - p1).cross(p3 - p1).norm();

                            double w = area;
                            ww += w;
                            u += w * avg;
                        }
                    }

                    if (ww > 0)
                    { // to avoid overflow (i.e., ww == 0)
                        u /= ww;
                        u -= p1;
                        // n = vnormal_[v];
                        // u -= n * dot(u, n);
                        updated_coords[n_id] = u;
                    }
                }
            }
        }

        // Update positions
        for (auto n_id : utility::range(surface_mesh.nodes().size()))
        {
            if (!is_node_on_corner[n_id])
            {
                simple_math::Vector3D new_point =
                  surface_mesh.nodes()[n_id].coord() + updated_coords[n_id];
                for (auto j : utility::range(3))
                {
                    surface_mesh.nodes()[n_id].point(j) = new_point(j);
                }
            }
        }

        if (GlobalOptions::get().write_debug_meshes)
        {
            mesh::write_vtk(surface_mesh,
                            GlobalOptions::get().output_dir + "/iteration_preproject_" +
                              std::to_string(iteration + 1) + ".vtk",
                            data);
        }

        // Project nodes back to the facets they should lie on
        for (auto n_id : utility::range(surface_mesh.nodes().size()))
        {
            if (!is_node_on_corner[n_id] && !is_node_on_seam[n_id])
            {
                simple_math::Vector3D new_point = surface_mesh.nodes()[n_id].coord();
                for (auto& facet_id : node_id_to_input_facet_ids_map[n_id])
                {
                    auto& igl_surface = igl_facets[facet_id];
                    new_point         = igl_closest_point(igl_surface, new_point);
                }


                for (auto j : utility::range(3))
                {
                    surface_mesh.nodes()[n_id].point(j) = new_point(j);
                }
            }
        }

        // Project nodes on seams back to the seams
        if (allow_smoothing_on_seams)
        {
            for (auto seam_id : utility::range(n_seams))
            {
                for (auto& n_id : seam_id_to_node_id_map[seam_id])
                {
                    if (!is_node_on_corner[n_id])
                    {
                        assert(node_to_neighbors_on_seam[n_id].size() == 2);
                        auto p = surface_mesh.nodes()[n_id].point();

                        auto& neighbor1 = node_to_neighbors_on_seam[n_id][0];
                        auto result1    = simple_math::closest_point_on_line_segment(
                          p,
                          original_positions[edges[neighbor1.edge_id][0]],
                          original_positions[edges[neighbor1.edge_id][1]]);

                        auto& neighbor2 = node_to_neighbors_on_seam[n_id][1];
                        auto result2    = simple_math::closest_point_on_line_segment(
                          p,
                          original_positions[edges[neighbor2.edge_id][0]],
                          original_positions[edges[neighbor2.edge_id][1]]);

                        Index curr_n_id    = n_id;
                        auto curr_neighbor = neighbor1;
                        auto curr_result   = result1;
                        if (result2.distance < result1.distance)
                        {
                            // The point is closer to the second neighbor so start iterating that
                            // direction instead of in the direction of the first neighbor.
                            curr_neighbor = neighbor2;
                            curr_result   = result2;
                        }

                        // Traverse the seam until an edge is reached that is farther away from the
                        // point.
                        while (true)
                        {
                            auto& next_neighbors = node_to_neighbors_on_seam[curr_neighbor.node_id];
                            if (next_neighbors.size() > 2)
                            {
                                // We traversed until we encounted a corner, so just use the last
                                // result
                                break;
                            }

                            auto next_neighbor = next_neighbors[0];
                            if (next_neighbors[0].node_id == curr_n_id)
                            {
                                next_neighbor = next_neighbors[1];
                            }
                            auto next_result = simple_math::closest_point_on_line_segment(
                              p,
                              original_positions[edges[next_neighbor.edge_id][0]],
                              original_positions[edges[next_neighbor.edge_id][1]]);

                            if (next_result.distance > curr_result.distance)
                            {
                                // The next segment is farther from the point than the last one, so
                                // stop traversing the seam
                                break;
                            }

                            curr_n_id     = curr_neighbor.node_id;
                            curr_neighbor = next_neighbor;
                            curr_result   = next_result;
                        }

                        surface_mesh.nodes()[n_id].point() = curr_result.closest_point;
                    }
                }
            }
        }


        // Calculate how much each node was nominally moved (before projections)
        simple_math::Vector3D avg_update = simple_math::Vector3D::Zero();
        for (auto n_id : utility::range(surface_mesh.nodes().size()))
        {
            for (auto j : utility::range(3))
            {
                avg_update(j) += (updated_coords[n_id][j] - avg_update(j)) / (double)(n_id + 1);
            }
        }
        std::cout << GlobalOptions::get().curr_status_padding << "   Vertices moved an average of "
                  << avg_update.norm() << std::endl;

        if (GlobalOptions::get().write_debug_meshes)
        {
            mesh::write_vtk(surface_mesh,
                            GlobalOptions::get().output_dir + "/iteration_" +
                              std::to_string(iteration + 1) + ".vtk",
                            data);
        }
    }


    if (GlobalOptions::get().write_debug_meshes)
    {
        mesh::write_vtk(surface_mesh,
                        GlobalOptions::get().output_dir + "/post_tangential_smoothing.vtk");
    }
}


std::vector<mesh::ConnIndex> calculate_neighbor_counts(
  const mesh::Mesh& surface_mesh,
  const mesh::ComponentIndex& component_to_smooth,
  const std::vector<bool>& is_node_constrained)
{
    std::vector<mesh::ConnIndex> neighbor_counts(surface_mesh.nodes().size());
    for (auto surf_element : surface_mesh.elements())
    {
        if (surf_element.component_id() == component_to_smooth)
        {
            for (auto node : surf_element.nodes())
            {
                if (!is_node_constrained[node.local_index()])
                {
                    auto& c = neighbor_counts[node.local_index()];
                    for (auto& node_2 : surf_element.nodes())
                    {
                        if (node != node_2)
                        {
                            c++;
                        }
                    }
                }
            }
        }
    }

    return neighbor_counts;
}


void optimize_tangential_smoothing(Triangulation& t,
                                   const SeamInfo& seams,
                                   const int n_iterations,
                                   const bool allow_smoothing_on_seams)
{
    using namespace simple_math::point_operators;

    std::cout << GlobalOptions::get().curr_status_padding
              << "* optimize_tangential_smoothing - Starting iterative optimization" << std::endl;

    if (GlobalOptions::get().write_debug_meshes)
    {
        mesh::write_vtk(to_mesh(t), GlobalOptions::get().output_dir + "/iteration_0.vtk");
        detail::write_seams(t, seams, GlobalOptions::get().output_dir + "/iteration_0_");
    }

    // Count # of facets and seams
    mesh::ComponentIndex n_facets =
      *std::max_element(t.tri_to_marker_map.begin(), t.tri_to_marker_map.end()) + 1;

    auto vert_to_tri_map = build_vert_to_tri_map(t);

    // Make a copy of the original node positions
    auto original_positions = t.vertices;

    // Create a map from each feature to the nodes that compose them
    std::vector<std::vector<mesh::LocalIndex>> facet_id_to_vert_id_map(n_facets);
    for (auto v_id : utility::indices_as<mesh::LocalIndex>(t.vertices))
    {
        for (auto tri_id : vert_to_tri_map.get_indices(v_id))
        {
            facet_id_to_vert_id_map[t.tri_to_marker_map[tri_id]].push_back(v_id);
        }
    }
    std::vector<bool> is_vert_on_seam(t.vertices.size(), false);
    std::vector<mesh::ComponentIndex> n_seam_elements_connected_to_vert(t.vertices.size(), 0);
    std::vector<std::vector<mesh::LocalIndex>> seam_id_to_vert_id_map(seams.n_seams);
    std::vector<std::vector<detail::AdjacentNode>> vert_to_neighbors_on_seam(t.vertices.size());
    for (auto feature_edge_id : seams.edges_on_seams)
    {
        auto seam_id = seams.edge_id_to_seam_id_map[feature_edge_id];

        auto v1_id             = t.edges[feature_edge_id][0];
        is_vert_on_seam[v1_id] = true;
        n_seam_elements_connected_to_vert[v1_id]++;
        seam_id_to_vert_id_map[seam_id].push_back(v1_id);

        auto v2_id             = t.edges[feature_edge_id][1];
        is_vert_on_seam[v2_id] = true;
        n_seam_elements_connected_to_vert[v2_id]++;
        seam_id_to_vert_id_map[seam_id].push_back(v2_id);

        vert_to_neighbors_on_seam[v1_id].push_back({ v2_id, feature_edge_id });
        vert_to_neighbors_on_seam[v2_id].push_back({ v1_id, feature_edge_id });
    }

    // Create a map from each node to the facets it should lie on
    std::vector<utility::VectorSet<unsigned int>> vert_id_to_input_facet_ids_map(t.vertices.size());
    for (auto tri_id : utility::indices(t.triangles))
    {
        for (auto v_id : t.triangles[tri_id])
        {
            vert_id_to_input_facet_ids_map[v_id].insert(t.tri_to_marker_map[tri_id]);
        }
    }

    // Figure out which nodes are on a corner based on how many seam elements are connected to
    // the node.  It is assumed that if it is greater than 2, then it a corner.
    std::vector<bool> is_vert_on_corner(t.vertices.size(), false);
    for (auto n_id : utility::indices(t.vertices))
    {
        if (n_seam_elements_connected_to_vert[n_id] > 2)
        {
            is_vert_on_corner[n_id] = true;
        }
    }

    // TODO build adjacency list for seam edges

    auto igl_facets = decompose_into_igl_facets(t);

    std::vector<uint32_t> is_vert_on_seam_data(is_vert_on_seam.size(), 0);
    for (auto i : utility::indices(is_vert_on_seam_data))
    {
        if (is_vert_on_seam[i])
        {
            is_vert_on_seam_data[i] = is_vert_on_seam[i];
        }
    }

    mesh::MeshAssociatedData data;
    data.add<mesh::DataFieldType::NODE, uint32_t>(
      mesh::DataField<uint32_t>::cref("is_on_seam", is_vert_on_seam_data));

    for (auto iteration : utility::range(n_iterations))
    {
        std::cout << GlobalOptions::get().curr_status_padding << " - Starting iteration "
                  << iteration << std::endl;

        std::vector<simple_math::Point> updated_coords(t.vertices.size());

        if (allow_smoothing_on_seams)
        {
            for (auto seam_id : utility::range(seams.n_seams))
            {
                for (auto& v_id : seam_id_to_vert_id_map[seam_id])
                {
                    if (!is_vert_on_corner[v_id])
                    {
                        simple_math::Vector3D u  = simple_math::Vector3D::Zero();
                        simple_math::Vector3D tv = simple_math::Vector3D::Zero();
                        double edge_length_sums  = 0;
                        int n_connected_edges    = 0;
                        auto p1                  = simple_math::point_as_vector(t.vertices[v_id]);

                        // Loop over points on feature connected to n_id via an edge
                        for (auto neighbor : vert_to_neighbors_on_seam[v_id])
                        {
                            auto p2 = simple_math::point_as_vector(t.vertices[neighbor.node_id]);

                            simple_math::Vector3D midpoint = 0.5 * (p1 + p2);
                            double distance                = (p2 - p1).norm();
                            edge_length_sums += distance;
                            u += distance * midpoint;

                            if (n_connected_edges == 0)
                            {
                                tv += (p2 - p1).normalized();
                                n_connected_edges++;
                            }
                            else
                            {
                                n_connected_edges++;
                                tv -= (p2 - p1).normalized();
                            }
                        }

                        assert(n_connected_edges == 2);

                        if (edge_length_sums > 0)
                        { // to avoid overflow (i.e., edge_length_sums == 0)
                            u /= edge_length_sums;
                            u -= p1;
                            tv                   = tv.normalized();
                            u                    = tv * u.dot(tv);
                            updated_coords[v_id] = simple_math::to_point(u);
                        }
                    }
                }
            }
        }

        for (auto facet_id : utility::range(n_facets))
        {
            for (auto& v_id : facet_id_to_vert_id_map[facet_id])
            {
                if (!is_vert_on_corner[v_id] && !is_vert_on_seam[v_id])
                {
                    simple_math::Vector3D u = simple_math::Vector3D::Zero();
                    double ww               = 0;
                    int c                   = 0;
                    auto p1                 = simple_math::point_as_vector(t.vertices[v_id]);

                    for (auto tri_id : vert_to_tri_map.get_indices(v_id))
                    {
                        if (t.tri_to_marker_map[tri_id] == facet_id)
                        {
                            auto p1 =
                              simple_math::point_as_vector(t.vertices[t.triangles[tri_id][0]]);
                            auto p2 =
                              simple_math::point_as_vector(t.vertices[t.triangles[tri_id][1]]);
                            auto p3 =
                              simple_math::point_as_vector(t.vertices[t.triangles[tri_id][2]]);
                            auto avg  = (p1 + p2 + p3) / 3.0;
                            auto area = (p2 - p1).cross(p3 - p1).norm();

                            double w = area;
                            ww += w;
                            u += w * avg;
                        }
                    }

                    if (ww > 0)
                    { // to avoid overflow (i.e., ww == 0)
                        u /= ww;
                        u -= p1;
                        // n = vnormal_[v];
                        // u -= n * dot(u, n);
                        updated_coords[v_id] = simple_math::to_point(u);
                    }
                }
            }
        }

        // Update positions
        for (auto v_id : utility::indices(t.vertices))
        {
            if (!is_vert_on_corner[v_id])
            {
                t.vertices[v_id] = t.vertices[v_id] + updated_coords[v_id];
            }
        }

        if (GlobalOptions::get().write_debug_meshes)
        {
            mesh::write_vtk(to_mesh(t),
                            GlobalOptions::get().output_dir + "/iteration_preproject_" +
                              std::to_string(iteration + 1) + ".vtk",
                            data);
        }

        // Project nodes back to the facets they should lie on
        for (auto v_id : utility::indices(t.vertices))
        {
            if (!is_vert_on_corner[v_id] && !is_vert_on_seam[v_id])
            {
                simple_math::Vector3D new_point = simple_math::point_as_vector(t.vertices[v_id]);
                for (auto& facet_id : vert_id_to_input_facet_ids_map[v_id])
                {
                    auto& igl_surface = igl_facets[facet_id];
                    new_point         = igl_closest_point(igl_surface, new_point);
                }

                t.vertices[v_id] = simple_math::to_point(new_point);
            }
        }

        // Project nodes on seams back to the seams
        if (allow_smoothing_on_seams)
        {
            for (auto seam_id : utility::range(seams.n_seams))
            {
                for (auto& v_id : seam_id_to_vert_id_map[seam_id])
                {
                    if (!is_vert_on_corner[v_id])
                    {
                        assert(vert_to_neighbors_on_seam[v_id].size() == 2);
                        auto p = t.vertices[v_id];

                        auto& neighbor1 = vert_to_neighbors_on_seam[v_id][0];
                        auto result1    = simple_math::closest_point_on_line_segment(
                          p,
                          original_positions[t.edges[neighbor1.edge_id][0]],
                          original_positions[t.edges[neighbor1.edge_id][1]]);

                        auto& neighbor2 = vert_to_neighbors_on_seam[v_id][1];
                        auto result2    = simple_math::closest_point_on_line_segment(
                          p,
                          original_positions[t.edges[neighbor2.edge_id][0]],
                          original_positions[t.edges[neighbor2.edge_id][1]]);

                        Index curr_n_id    = v_id;
                        auto curr_neighbor = neighbor1;
                        auto curr_result   = result1;
                        if (result2.distance < result1.distance)
                        {
                            // The point is closer to the second neighbor so start iterating that
                            // direction instead of in the direction of the first neighbor.
                            curr_neighbor = neighbor2;
                            curr_result   = result2;
                        }

                        // Traverse the seam until an edge is reached that is farther away from the
                        // point.
                        while (true)
                        {
                            auto& next_neighbors = vert_to_neighbors_on_seam[curr_neighbor.node_id];
                            if (next_neighbors.size() > 2)
                            {
                                // We traversed until we encounted a corner, so just use the last
                                // result
                                break;
                            }

                            auto next_neighbor = next_neighbors[0];
                            if (next_neighbors[0].node_id == curr_n_id)
                            {
                                next_neighbor = next_neighbors[1];
                            }
                            auto next_result = simple_math::closest_point_on_line_segment(
                              p,
                              original_positions[t.edges[next_neighbor.edge_id][0]],
                              original_positions[t.edges[next_neighbor.edge_id][1]]);

                            if (next_result.distance > curr_result.distance)
                            {
                                // The next segment is farther from the point than the last one, so
                                // stop traversing the seam
                                break;
                            }

                            curr_n_id     = curr_neighbor.node_id;
                            curr_neighbor = next_neighbor;
                            curr_result   = next_result;
                        }

                        t.vertices[v_id] = curr_result.closest_point;
                    }
                }
            }
        }


        // Calculate how much each node was nominally moved (before projections)
        simple_math::Vector3D avg_update = simple_math::Vector3D::Zero();
        for (auto n_id : utility::range(t.vertices.size()))
        {
            for (auto j : utility::range(3))
            {
                avg_update(j) += (updated_coords[n_id][j] - avg_update(j)) / (double)(n_id + 1);
            }
        }
        std::cout << GlobalOptions::get().curr_status_padding << "   Vertices moved an average of "
                  << avg_update.norm() << std::endl;

        if (GlobalOptions::get().write_debug_meshes)
        {
            mesh::write_vtk(to_mesh(t),
                            GlobalOptions::get().output_dir + "/iteration_" +
                              std::to_string(iteration + 1) + ".vtk",
                            data);
        }
    }


    if (GlobalOptions::get().write_debug_meshes)
    {
        mesh::write_vtk(to_mesh(t),
                        GlobalOptions::get().output_dir + "/post_tangential_smoothing.vtk");
    }
}

namespace detail
{

std::vector<mesh::ConnIndex> calculate_neighbor_counts(
  const mesh::Mesh& surface_mesh,
  const mesh::ComponentIndex& component_to_smooth,
  const std::vector<bool>& is_node_constrained)
{
    std::vector<mesh::ConnIndex> neighbor_counts(surface_mesh.nodes().size());
    for (auto surf_element : surface_mesh.elements())
    {
        if (surf_element.component_id() == component_to_smooth)
        {
            for (auto node : surf_element.nodes())
            {
                if (!is_node_constrained[node.local_index()])
                {
                    auto& c = neighbor_counts[node.local_index()];
                    for (auto& node_2 : surf_element.nodes())
                    {
                        if (node != node_2)
                        {
                            c++;
                        }
                    }
                }
            }
        }
    }

    return neighbor_counts;
}

/// @see https://graphics.stanford.edu/courses/cs468-01-fall/Papers/taubin-smoothing.pdf
void update_delta(const mesh::Mesh& surface_mesh,
                  const mesh::ComponentIndex& component_to_smooth,
                  const std::vector<bool>& is_node_constrained,
                  const std::vector<mesh::ConnIndex> neighbor_counts,
                  std::vector<simple_math::Vector3D>& delta_v)
{
    for (auto& d : delta_v)
    {
        d.setZero();
    }

    for (auto surf_element : surface_mesh.elements())
    {
        if (surf_element.component_id() == component_to_smooth)
        {
            for (auto& node : surf_element.nodes())
            {
                if (!is_node_constrained[node.local_index()])
                {
                    auto& d = delta_v[node.local_index()];
                    for (auto& node_2 : surf_element.nodes())
                    {
                        if (node != node_2)
                        {
                            d += node_2.coord() - node.coord();
                        }
                    }
                }
            }
        }
    }
    for (auto i : utility::range(delta_v.size()))
    {
        if (neighbor_counts[i] > 0)
        {
            delta_v[i] /= neighbor_counts[i];
        }
    }
}
}

/// @see https://graphics.stanford.edu/courses/cs468-01-fall/Papers/taubin-smoothing.pdf
void optimize_taubin_smoothing(mesh::Mesh& surface_mesh,
                               const mesh::ComponentIndex& component_to_smooth,
                               const std::vector<bool>& is_node_constrained,
                               const double lambda,
                               const double mu,
                               const int n_smoothing_iterations)
{
    auto neighbor_counts =
      calculate_neighbor_counts(surface_mesh, component_to_smooth, is_node_constrained);
    std::vector<simple_math::Vector3D> delta_v(surface_mesh.nodes().size());
    for (auto n : utility::range(n_smoothing_iterations))
    {
        detail::update_delta(surface_mesh,
                             component_to_smooth,
                             is_node_constrained,
                             neighbor_counts,
                             delta_v);
        for (auto i : utility::indices_as<Index>(delta_v))
        {
            if (!is_node_constrained[i])
            {
                auto& node = surface_mesh.nodes()[i];
                for (auto j : utility::range(3))
                {
                    node.point(j) += lambda * delta_v[i][j];
                }
            }
        }

        detail::update_delta(surface_mesh,
                             component_to_smooth,
                             is_node_constrained,
                             neighbor_counts,
                             delta_v);
        for (auto i : utility::range(delta_v.size()))
        {
            if (!is_node_constrained[i])
            {
                auto& node = surface_mesh.nodes()[i];
                for (auto j : utility::range(3))
                {
                    node.point(j) += mu * delta_v[i][j];
                }
            }
        }
    }
}

namespace detail
{
int find_external_vertex_position(const Triangulation& t, const Index& tri_id, const Index& edge_id)
{
    auto& edge = t.edges[edge_id];
    for (auto i : utility::range(3))
    {
        const auto& v_id = t.triangles[tri_id][i];
        if (v_id != edge[0] && v_id != edge[1])
        {
            return i;
        }
    }
    throw _UTIL_MSG_EXCEPT(
      "Bad triangle and edge combination given. The triangle may be collapsed.");
}

int find_external_vertex_valence(const Triangulation& t,
                                 const Index& tri_id,
                                 const Index& edge_id,
                                 const std::vector<int>& vertex_valences)
{
    auto i = find_external_vertex_position(t, tri_id, edge_id);
    return vertex_valences[t.triangles[tri_id][i]];
}

// Returns 0 if aligned and 1 if not aligned if t1 and the opposite if t2 (not t1).
int get_aligned_case(const Triangulation& t,
                     const bool& is_t1,
                     const Index& tri_id,
                     const Index& edge_id,
                     const int& edge_position)
{
    auto edge_info = mesh::get_edge_info(mesh::ElementType::T2D3N);
    if (t.edges[edge_id][0] == t.triangles[tri_id][edge_info[edge_position].connectivity[0]])
    {
        if (is_t1)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else if (t.edges[edge_id][1] == t.triangles[tri_id][edge_info[edge_position].connectivity[0]])
    {
        if (is_t1)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        throw _UTIL_MSG_EXCEPT(
          "Bad triangle and edge combination given. The edge may not be in the triangle.");
    }
}

void flip_edge(Triangulation& t, const Index& edge_to_flip_id)
{
    auto& connected_tris = t.edge_to_tris_map[edge_to_flip_id];
    assert(connected_tris.size() == 2);
    assert(t.tri_to_marker_map[connected_tris[0].first] ==
           t.tri_to_marker_map[connected_tris[1].first]);

    // We do the following transform to flip e = { 0, 2 }  ->  { 1, 3 }
    //
    //   3 e2 2         3 e2 2
    //   o----o         o----o
    //   |t1 /|         |\ t2|
    //   |  / |         | \e |
    // e3| /e |e1 ==> e3|  \ |e1
    //   |/ t2|         |t1 \|
    //   o----o         o----o
    //   0 e0 1         0 e0 1
    // Note: for the transform, we will make the following assumptions:
    //      1. Keep the position of the flipped edge in each triangle constant
    //      2. t1 will keep the original start point of the edge to flip
    //      3. t2 will keep the original end point of the edge to flip

    // Maps for t1 transformation
    //
    // Case edge is position 0:
    //           Aligned w/ Edge    Not Aligned w/ Edge
    //   t1 conn   { 0, 2, 3 }         { 2, 0, 3 }
    //   t1 edges { e, e2, e3 }       { e, e3, e2 }
    //
    // Case edge is position 1:
    //           Aligned w/ Edge    Not Aligned w/ Edge
    //   t1 conn   { 3, 0, 2 }         { 3, 2, 0 }
    //   t1 edges { e3, e, e2 }       { e2, e, e3 }
    //
    // Case edge is position 2:
    //           Aligned w/ Edge    Not Aligned w/ Edge
    //   t1 conn   { 2, 3, 0 }         { 0, 3, 2 }
    //   t1 edges { e2, e3, e }       { e3, e2, e }
    //
    // Note: maps for t2 transformation is same as t1, just reverse which is used for
    // aligned/unaligned cases. get_aligned_case adjusts for the reversal.
    //
    // Case edge is position 0:
    //         Not Aligned w/ Edge    Aligned w/ Edge
    //   t2 conn   { 2, 0, 1 }          { 0, 2, 1 }
    //   t2 edges { e, e0, e1 }        { e, e1, e0 }
    //
    // Case edge is position 1:
    //         Not Aligned w/ Edge    Aligned w/ Edge
    //   t2 conn   { 1, 2, 0 }          { 1, 0, 2 }
    //   t2 edges { e1, e, e0 }        { e0, e, e1 }
    //
    // Case edge is position 2:
    //         Not Aligned w/ Edge    Aligned w/ Edge
    //   t2 conn   { 0, 1, 2 }         { 2, 1, 0 }
    //   t2 edges { e0, e1, e }       { e1, e0, e }

    // The location of the external vertex in the triangle's original connectivity (used for v3 and
    // v1) Note: whether the triangle and edge align does not affect this map Remember: the index is
    // the "case" which is the position of the edge to flip in the triangle
    static int ext_v_pos[3] = {
        2, // case 0
        0, // case 1
        1  // case 2
    };
    // Indicates the position of e2 in t1 or e0 in t2 (before the flip)
    static int ext_e_pos[2][3] = { {
                                     // t1 is aligned with edge
                                     1, // case 0
                                     2, // case 1
                                     0  // case 2
                                   },
                                   {
                                     // t1 is not aligned with edge
                                     2, // case 0
                                     0, // case 1
                                     1  // case 2
                                   } };
    // Indicates the position of e3 in t1 or e1 in t2 (before the flip)
    static int reused_e_pos[2][3] = { {
                                        // t1 is aligned with edge
                                        2, // case 0
                                        0, // case 1
                                        1  // case 2
                                      },
                                      {
                                        // t2 is not aligned with edge
                                        1, // case 0
                                        2, // case 1
                                        0  // case 2
                                      } };

    auto& t1_id = connected_tris[0].first;
    auto& t2_id = connected_tris[1].first;

    const auto& t1_case = connected_tris[0].second;
    auto t1_align_i     = get_aligned_case(t, true, t1_id, edge_to_flip_id, t1_case);
    const auto& t2_case = connected_tris[1].second;
    auto t2_align_i     = get_aligned_case(t, false, t2_id, edge_to_flip_id, t2_case);

    auto v0_id = t.edges[edge_to_flip_id][0];
    auto v1_id = t.triangles[t2_id][ext_v_pos[t2_case]];
    auto v2_id = t.edges[edge_to_flip_id][1];
    auto v3_id = t.triangles[t1_id][ext_v_pos[t1_case]];

    auto e0_id = t.tri_to_edges_map[t2_id][ext_e_pos[t2_align_i][t2_case]];
    auto e1_id = t.tri_to_edges_map[t2_id][reused_e_pos[t2_align_i][t2_case]];
    auto e2_id = t.tri_to_edges_map[t1_id][ext_e_pos[t1_align_i][t1_case]];
    auto e3_id = t.tri_to_edges_map[t1_id][reused_e_pos[t1_align_i][t1_case]];

    auto replace_edge = [&t](const Index& edge_id,
                             const Index& tri_to_find,
                             const Index& tri_id,
                             const Index& edge_position)
    {
        for (auto& tri : t.edge_to_tris_map[edge_id])
        {
            if (tri.first == tri_to_find)
            {
                tri.first  = tri_id;
                tri.second = edge_position;
                return;
            }
        }
        throw _UTIL_MSG_EXCEPT("Triangle " + std::to_string(tri_to_find) +
                               " did not appear in map for edge " + std::to_string(edge_id));
    };

    t.edges[edge_to_flip_id] = { v1_id, v3_id };

    // Update for tri 1
    if (t1_align_i == 0) // aligned
    {
        t.triangles[t1_id]        = { v1_id, v3_id, v0_id };
        t.tri_to_edges_map[t1_id] = { edge_to_flip_id, e3_id, e0_id };
        replace_edge(edge_to_flip_id, t1_id, t1_id, 0);
        replace_edge(e3_id, t1_id, t1_id, 1);
        replace_edge(e0_id, t2_id, t1_id, 2);
    }
    else // not aligned
    {
        t.triangles[t1_id]        = { v3_id, v1_id, v2_id };
        t.tri_to_edges_map[t1_id] = { edge_to_flip_id, e0_id, e3_id };
        replace_edge(edge_to_flip_id, t1_id, t1_id, 0);
        replace_edge(e0_id, t2_id, t1_id, 1);
        replace_edge(e3_id, t1_id, t1_id, 2);
    }
    // Update for tri 2
    if (t2_align_i == 0) // not aligned, pointing out of plane of figure
    {
        t.triangles[t2_id]        = { v3_id, v1_id, v2_id };
        t.tri_to_edges_map[t2_id] = { edge_to_flip_id, e1_id, e2_id };
        replace_edge(edge_to_flip_id, t2_id, t2_id, 0);
        replace_edge(e1_id, t2_id, t2_id, 1);
        replace_edge(e2_id, t1_id, t2_id, 2);
    }
    else // aligned, pointing into plane of figure
    {
        t.triangles[t2_id]        = { v1_id, v3_id, v2_id };
        t.tri_to_edges_map[t2_id] = { edge_to_flip_id, e2_id, e1_id };
        replace_edge(edge_to_flip_id, t2_id, t2_id, 0);
        replace_edge(e2_id, t1_id, t2_id, 1);
        replace_edge(e1_id, t2_id, t2_id, 2);
    }
}
}

void optimize_edge_flips(Triangulation& t, const int n_passes)
{
    check_tri_edge_relations(t);

    // Calculate the number of edges attached to each node
    std::vector<int> vertex_valences(t.vertices.size(), 0);
    std::vector<utility::DynamicArray<Index, mesh::ConnIndex>> node_to_edge_map(t.vertices.size());
    for (auto edge_id : utility::indices(t.edges))
    {
        auto& edge = t.edges[edge_id];
        vertex_valences[edge[0]]++;
        vertex_valences[edge[1]]++;
        node_to_edge_map[edge[0]].push_back(edge_id);
    }

    std::cout << GlobalOptions::get().curr_status_padding
              << "* optimize_edge_flips - Starting edge flips, current max valence is "
              << *std::max_element(vertex_valences.begin(), vertex_valences.end()) << std::endl;

    for (auto pass : utility::range(n_passes))
    {
        // Find the permutation that sorts the nodes from highest to lowest valence.
        // In other words, sorted_permutation[0] is the node ID with highest.
        auto sorted_permutation = utility::find_sorted_permutation(vertex_valences);
        std::reverse(sorted_permutation.begin(), sorted_permutation.end());

        auto pass_start_max_valence = vertex_valences[sorted_permutation[0]];

        // Fig 1. Example that will be used for comments in this function
        //
        //   nodes      edges    triangles
        // 6---7---8  o---o---o  o---o---o
        // |\  |  /|  |\  |  /|  |\ 5|6 /|
        // | \ | / |  | 5 6 7 |  | \ | / |
        // |  \|/  |  |  \|/  |  |4 \|/ 7|
        // 5---0---1  o-4-o-0-o  o---o---o
        // |  /|\  |  |  /|\  |  |3 /|\ 0|
        // | / | \ |  | 3 2 1 |  | / | \ |
        // |/  |  \|  |/  |  \|  |/ 2|1 \|
        // 4---3---2  o---o---o  o---o---o
        //

        // Iterate until all vertices with valence higher than 7 have been visited
        std::size_t curr_i = 0;
        while (vertex_valences[sorted_permutation[curr_i]] > 7)
        {
            auto v_id       = sorted_permutation[curr_i];
            auto& v_valence = vertex_valences[v_id];

            // Let "external nodes" for a given edge be the corner nodes of the triangles connected
            // to the edge that are not in the edge.  For Fig 1, if v_id=0, then 1 and 7 are the
            // external nodes for edge 7.
            //
            // Find the max valence of the external nodes for each edge, since we want to flip edge
            // that have a low "max external valence".
            auto orig_connected_edges = node_to_edge_map[v_id];
            std::vector<int> initial_max_external_valences(v_valence, 1000);
            std::vector<std::array<Index, 2>> external_vertices(v_valence, { 0, 0 });
            std::vector<bool> is_connected_edge_locked(v_valence, false);
            for (auto i : utility::range(orig_connected_edges.size()))
            {
                auto& edge_id = orig_connected_edges[i];
                auto& edge    = t.edges[edge_id];

                // Do not modify an edge that is connected to more than 2 triangles or triangles
                // with different markers (since we would flip an edge on a seam).
                if (t.edge_to_tris_map[edge_id].size() != 2 ||
                    t.tri_to_marker_map[t.edge_to_tris_map[edge_id][0].first] !=
                      t.tri_to_marker_map[t.edge_to_tris_map[edge_id][1].first])
                {
                    is_connected_edge_locked[i] = true;
                    continue;
                }

                auto& t1_id = t.edge_to_tris_map[edge_id][0].first;
                auto i_1    = detail::find_external_vertex_position(t, t1_id, edge_id);

                auto& t2_id = t.edge_to_tris_map[edge_id][1].first;
                auto i_2    = detail::find_external_vertex_position(t, t2_id, edge_id);

                external_vertices[i] = { t.triangles[t1_id][i_1], t.triangles[t2_id][i_2] };

                initial_max_external_valences[i] =
                  std::min(vertex_valences[external_vertices[i][0]],
                           vertex_valences[external_vertices[i][1]]);
            }

            // Permutation that sorts the neighbors from those with the lowest "max external
            // valence" to the highest. Note: their valences might change as edges are flipped (this
            // is the permutation based on the initial values before any flips).
            auto neighbor_sorted_p =
              utility::find_sorted_permutation(initial_max_external_valences);
            for (auto i : utility::indices_as<Index>(orig_connected_edges))
            {
                auto& sorted_i = neighbor_sorted_p[i];
                auto edge_id   = orig_connected_edges[sorted_i];

                // If the current vertex has a higher valence than both the external nodes of the
                // currently considered edge and the edge is not locked, then flip the edge.
                // Note: use the up-to-date valences.
                if (v_valence > vertex_valences[external_vertices[sorted_i][0]] + 1 &&
                    v_valence > vertex_valences[external_vertices[sorted_i][1]] + 1 &&
                    !is_connected_edge_locked[sorted_i])
                {
                    auto orig_edge = t.edges[edge_id];
                    detail::flip_edge(t, edge_id); // changes connectivity of edge

                    // Update valences
                    vertex_valences[orig_edge[0]]--;
                    vertex_valences[orig_edge[1]]--;
                    vertex_valences[t.edges[edge_id][0]]++;
                    vertex_valences[t.edges[edge_id][1]]++;

                    // Update node-to-edge map
                    node_to_edge_map[orig_edge[0]].remove_value(edge_id);
                    node_to_edge_map[orig_edge[1]].remove_value(edge_id);
                    node_to_edge_map[t.edges[edge_id][0]].push_back(edge_id);
                    node_to_edge_map[t.edges[edge_id][1]].push_back(edge_id);
                }
            }

            // auto max_iter = std::max_element(vertex_valences.begin(), vertex_valences.end());
            // std::cout << " - Flipped edges for v " << v_id << " max valence is " << *max_iter
            //           << " (at v " << std::distance(vertex_valences.begin(), max_iter) << ")"
            //           << std::endl;

            curr_i++;
        }

        auto pass_end_max_valence =
          *std::max_element(vertex_valences.begin(), vertex_valences.end());
        auto pass_valence_reduction = pass_start_max_valence - pass_end_max_valence;
        std::cout << GlobalOptions::get().curr_status_padding << " - Pass " << pass
                  << " reduced max valence by " << pass_valence_reduction << " to "
                  << pass_end_max_valence << std::endl;
        if (pass_valence_reduction == 0)
        {
            std::cout << GlobalOptions::get().curr_status_padding
                      << " - Stopping since valence was not reduced in last pass." << std::endl;
            break;
        }
    }
}

}