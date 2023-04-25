#include "comp_mesh_core/Optimize.h"

#include "comp_mesh_core/GlobalOptions.h"
#include "comp_mesh_core/IGLHelpers.h"
#include "mesh/algorithms/IO.h"
#include "mesh/relations/NodeToElementMap.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"
#include "utility/VectorSet.h"

namespace composite_mesher
{
namespace detail
{

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
}

struct AdjacentNode
{
    Index node_id;
    Index edge_id;
};

void optimize_tangential_smoothing(mesh::Mesh& surface_mesh,
                                   const std::vector<std::array<mesh::LocalIndex, 2>>& edges,
                                   const mesh::LocalIndex& n_seams,
                                   const std::vector<mesh::LocalIndex>& edges_on_seams,
                                   const std::vector<mesh::LocalIndex>& edge_id_to_seam_id_map,
                                   const bool allow_smoothing_on_seams)
{
    std::cout << "* optimize_tangential_smoothing - Starting iterative optimization" << std::endl;

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
    std::vector<std::vector<AdjacentNode>> node_to_neighbors_on_seam(surface_mesh.nodes().size());
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
        std::cout << "  - Starting iteration " << iteration << std::endl;

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

        // TODO project nodes on seams back to the seams
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
        update_delta(surface_mesh,
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
                    node.point(j) += lambda * delta_v[i][j];
                }
            }
        }

        update_delta(surface_mesh,
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


#ifdef adfadfasdf

#pragma region approach_1

// Note: this suffers from a lack of robustness
void optimize_refinement(SurfaceMeshWithFeatures& feature_surface, const double& target_edge_length)
{
    mesh::write_vtk(feature_surface.surface_mesh, "optimize_refinement_input.vtk");

    mesh::write_vtk(get_seam_mesh(feature_surface.seams, feature_surface.surface_mesh.nodes()),
                    GlobalOptions::get().output_dir + "/optimize_input_seams.vtk");

    ////////////////////////////////////////////////////////////////////////////////
    // Convert triangular mesh to polygons and refine edges

    // Swap element types for polygons
    for (auto e : feature_surface.surface_mesh.elements())
    {
        e.type() = mesh::ElementType::POLYGON;
    }

    // TODO fix max length
    auto orig_edge_hash_to_new_conn_map =
      mesh::refine_all_polygon_segments(feature_surface.surface_mesh, target_edge_length);

    // Update the seam info for the refined segments
    const auto orig_seam_edge_size = feature_surface.seams.edge_id_to_connectivity_map.size();
    for (auto i : utility::range(orig_seam_edge_size))
    {
        auto& orig_conn = feature_surface.seams.edge_id_to_connectivity_map[i];
        auto hash       = utility::hash_sorted(orig_conn);
        auto find_iter  = orig_edge_hash_to_new_conn_map.find(hash);
        if (find_iter != orig_edge_hash_to_new_conn_map.end())
        {
            auto new_conn = find_iter->second; // copy intentional
            // They are ordered the opposite direction, so flip it
            if (orig_conn.first != new_conn[0])
            {
                std::reverse(new_conn.begin(), new_conn.end());
            }
            orig_conn.second   = new_conn[1];
            const auto seam_id = feature_surface.seams.edge_id_to_seam_id_map[i];
            // Note: orig_conn falls out of scope if this loop is entered
            for (auto j : utility::range(std::size_t(1), new_conn.size() - 1))
            {
                feature_surface.seams.edge_id_to_connectivity_map.push_back(
                  { new_conn[j], new_conn[j + 1] });
                feature_surface.seams.edge_id_to_seam_id_map.push_back(seam_id);
            }
        }
    }

    mesh::write_vtk(feature_surface.surface_mesh, "optimize_refinement_refined_edges.vtk");

    mesh::write_vtk(get_seam_mesh(feature_surface.seams, feature_surface.surface_mesh.nodes()),
                    GlobalOptions::get().output_dir + "/optimize_0_seams.vtk");

    ////////////////////////////////////////////////////////////////////////////////
    // Collapse short edges and remove degenerate polygons
    /*
    auto collapse_result =
      mesh::collapse_short_polygon_segments(feature_surface.surface_mesh,
                                            0.5 * target_edge_length,
                                            feature_surface.seams.edge_id_to_connectivity_map,
                                            feature_surface.seams.edge_id_to_seam_id_map);

    // Update interface information for the elements removed from collapsed segments that result
    // in two or less points.
    std::sort(collapse_result.elements_removed.begin(), collapse_result.elements_removed.end());
    std::reverse(collapse_result.elements_removed.begin(),
    collapse_result.elements_removed.end()); for (auto& i : collapse_result.elements_removed)
    {
        feature_surface.triangle_id_to_input_facet_ids_map.erase(
          feature_surface.triangle_id_to_input_facet_ids_map.begin() + i);
    }
    // Update seam information for the collapsed segments.
    std::unordered_map<std::size_t, mesh::LocalIndex> edge_hash_to_edge_id_map;
    for (auto i : utility::indices(feature_surface.seams.edge_id_to_connectivity_map))
    {
        auto& conn                     = feature_surface.seams.edge_id_to_connectivity_map[i];
        auto hash                      = utility::hash_sorted(conn);
        edge_hash_to_edge_id_map[hash] = i;
    };
    std::vector<mesh::LocalIndex> edge_ids_to_remove;
    for (auto& p : collapse_result.node_combinations)
    {
        auto hash = utility::hash_sorted(p);
        edge_ids_to_remove.push_back(edge_hash_to_edge_id_map[hash]);
    }
    std::sort(edge_ids_to_remove.begin(), edge_ids_to_remove.end());
    std::reverse(edge_ids_to_remove.begin(), edge_ids_to_remove.end());
    for (auto& i : edge_ids_to_remove)
    {
        feature_surface.seams.edge_id_to_seam_id_map.erase(
          feature_surface.seams.edge_id_to_seam_id_map.begin() + i);
        feature_surface.seams.edge_id_to_connectivity_map.erase(
          feature_surface.seams.edge_id_to_connectivity_map.begin() + i);
    }

    mesh::write_vtk(feature_surface.surface_mesh, "optimize_refinement_collapsed_edges.vtk");

    mesh::write_vtk(get_seam_mesh(feature_surface.seams, feature_surface.surface_mesh.nodes()),
                    GlobalOptions::get().output_dir + "/optimize_1_seams.vtk");
    */
    ////////////////////////////////////////////////////////////////////////////////
    // Create quality triangulation from the polygonal surface

    const double max_area = target_edge_length * target_edge_length;
    auto orig_n_elements  = feature_surface.surface_mesh.elements().size();

    // Loop over elements in surface mesh and submesh polygons with more than 4 nodes
    std::vector<mesh::LocalIndex> elements_to_remove;
    for (auto polygon_element_index :
         utility::range(feature_surface.surface_mesh.elements().size()))
    {
        auto& element = feature_surface.surface_mesh.elements()[polygon_element_index];
        if (element.connectivity().size() > 3)
        {
            std::vector<std::vector<mesh::LocalIndex>> boundary_segments;
            auto& conn = element.connectivity();
            for (auto segment : mesh::polygon_segments(element.connectivity()))
            {
                boundary_segments.push_back({ conn[segment.first], conn[segment.second] });
            }
            auto new_triangles =
              triangulate_polygon(feature_surface.surface_mesh, element, max_area, false);
            elements_to_remove.push_back(element.local_index());
            auto facet_ids =
              feature_surface.triangle_id_to_input_facet_ids_map[element.local_index()];
            for (auto& t : new_triangles)
            {
                feature_surface.triangle_id_to_input_facet_ids_map.push_back(facet_ids);
            }
        }
        else if (element.connectivity().size() == 3)
        {
            // Convert polygons to trianlges
            element.type() = mesh::ElementType::T2D3N;
        }
        else
        {
            throw _UTIL_MSG_EXCEPT("Bad element in surface.");
        }
    }

    auto n_new_triangles = feature_surface.surface_mesh.elements().size() - orig_n_elements;

    // Note: elements_to_remove is naturally sorted
    for (auto i : utility::range(elements_to_remove.size()))
    {
        feature_surface.surface_mesh.elements().remove_stable(elements_to_remove[i] - i);
        feature_surface.triangle_id_to_input_facet_ids_map.erase(
          feature_surface.triangle_id_to_input_facet_ids_map.begin() + elements_to_remove[i] - i);
    }

    mesh::write_vtk(feature_surface.surface_mesh, "optimize_refinement_triangulated.vtk");

    mesh::write_vtk(get_seam_mesh(feature_surface.seams, feature_surface.surface_mesh.nodes()),
                    GlobalOptions::get().output_dir + "/optimize_2_seams.vtk");
}

static void optimize_corefinement_step(SurfaceMeshWithFeatures& feature_surface)
{
    using namespace simple_math::point_operators;

    mesh::Mesh new_surface_mesh;
    for (auto n : feature_surface.surface_mesh.nodes())
    {
        new_surface_mesh.nodes().add(n.point());
    }
    std::vector<utility::VectorSet<unsigned int>> new_triangle_id_to_input_facet_ids_map;

    std::unordered_map<std::size_t, mesh::LocalIndex> orig_edge_hash_to_new_midpoint_map;
    for (auto e : feature_surface.surface_mesh.elements())
    {
        std::vector<mesh::LocalIndex> polygon_conn;

        auto conn = e.connectivity();
        for (auto s : mesh::polygon_segments(conn))
        {
            auto& n1_id = conn[s.first];
            auto& n2_id = conn[s.second];
            polygon_conn.push_back(n1_id);

            auto hash      = utility::hash_sorted(n1_id, n2_id);
            auto find_iter = orig_edge_hash_to_new_midpoint_map.find(hash);
            if (find_iter == orig_edge_hash_to_new_midpoint_map.end())
            {
                // This edge hasn't been visited, so create a new midpoint for this edge.
                auto midpoint = (feature_surface.surface_mesh.nodes()[n1_id].coord() +
                                 feature_surface.surface_mesh.nodes()[n2_id].coord()) /
                                2.0;

                new_surface_mesh.nodes().add(simple_math::to_point(midpoint));
                auto midpoint_id = new_surface_mesh.nodes().size() - 1;
                polygon_conn.push_back(midpoint_id);
                orig_edge_hash_to_new_midpoint_map[hash] = midpoint_id;
            }
            else
            {
                polygon_conn.push_back(find_iter->second);
            }
        }

        //     x 4
        //    /3\
        // 5 o---o 3
        //  /0\1/2\
        // x---o---x
        // 0   1   2
        new_surface_mesh.elements().add(mesh::ElementType::T2D3N,
                                        { polygon_conn[0], polygon_conn[1], polygon_conn[5] });
        new_surface_mesh.elements().add(mesh::ElementType::T2D3N,
                                        { polygon_conn[1], polygon_conn[3], polygon_conn[5] });
        new_surface_mesh.elements().add(mesh::ElementType::T2D3N,
                                        { polygon_conn[1], polygon_conn[2], polygon_conn[3] });
        new_surface_mesh.elements().add(mesh::ElementType::T2D3N,
                                        { polygon_conn[3], polygon_conn[4], polygon_conn[5] });
        auto& facet_ids = feature_surface.triangle_id_to_input_facet_ids_map[e.local_index()];
        for (auto i : utility::range(4))
        {
            new_triangle_id_to_input_facet_ids_map.push_back(facet_ids);
        }
    }

    // Update the seam info for the refined edges

    std::unordered_map<std::size_t, mesh::LocalIndex> feature_edge_hash_to_id_map;
    for (auto i : utility::indices(feature_surface.seams.edge_id_to_connectivity_map))
    {
        auto& conn                        = feature_surface.seams.edge_id_to_connectivity_map[i];
        auto hash                         = utility::hash_sorted(conn);
        feature_edge_hash_to_id_map[hash] = i;
    };

    auto& feature_edges         = feature_surface.seams.edge_id_to_connectivity_map;
    auto& edge_id_to_feature_id = feature_surface.seams.edge_id_to_seam_id_map;
    for (auto& v : orig_edge_hash_to_new_midpoint_map)
    {
        auto feature_edge_id                  = feature_edge_hash_to_id_map[v.first];
        auto tmp                              = feature_edges[feature_edge_id].second;
        feature_edges[feature_edge_id].second = v.second;
        feature_edges.push_back({ v.second, tmp });
        edge_id_to_feature_id.push_back(edge_id_to_feature_id[feature_edge_id]);
    }

    feature_surface.surface_mesh = std::move(new_surface_mesh);
    feature_surface.triangle_id_to_input_facet_ids_map =
      std::move(new_triangle_id_to_input_facet_ids_map);
}

void optimize_corefinement(SurfaceMeshWithFeatures& feature_surface,
                           const double target_mean_edge_length)
{
    auto stats = mesh::calc_segment_length_stats(feature_surface.surface_mesh);
    while (stats.mean > target_mean_edge_length)
    {
        optimize_corefinement_step(feature_surface);
        stats = mesh::calc_segment_length_stats(feature_surface.surface_mesh);
    }
}

mesh::Mesh tetrahedralize(SurfaceMeshWithFeatures& feature_surface)
{
    mesh::Mesh surface_mesh;
    surface_mesh.copy(feature_surface.surface_mesh);

    // Remove all elements but triangular elements
    std::vector<mesh::LocalIndex> elements_to_remove;
    for (auto& e : surface_mesh.elements())
    {
        if (e.type() != mesh::ElementType::T2D3N)
        {
            elements_to_remove.push_back(e.local_index());
        }
    }
    surface_mesh.elements().remove_stable(elements_to_remove);

    auto volume_mesh = mesh::tetrahedralize_surface(surface_mesh);

    // Assign materials via winding number?

    return volume_mesh;
}

#pragma endregion

#pragma region approach_2

static void add_directed_edge(std::unordered_map<mesh::LocalIndex, std::array<mesh::LocalIndex, 2>>&
                                boundary_node_adjacency_list,
                              const mesh::LocalIndex& start,
                              const mesh::LocalIndex& end)
{
    // Add the outward edge from start
    auto find_iter = boundary_node_adjacency_list.find(start);
    if (find_iter == boundary_node_adjacency_list.end())
    {
        boundary_node_adjacency_list[start] = { end, std::numeric_limits<mesh::LocalIndex>::max() };
    }
    else
    {
        assert(find_iter->second[0] == std::numeric_limits<mesh::LocalIndex>::max());
        find_iter->second[0] = end;
    }

    // Add the inward edge into end
    find_iter = boundary_node_adjacency_list.find(end);
    if (find_iter == boundary_node_adjacency_list.end())
    {
        boundary_node_adjacency_list[end] = { std::numeric_limits<mesh::LocalIndex>::max(), start };
    }
    else
    {
        assert(find_iter->second[1] == std::numeric_limits<mesh::LocalIndex>::max());
        find_iter->second[1] = start;
    }
}

mesh::Mesh tetrahedralize_via_direct_tetgen(const std::vector<bool>& is_facet_planar,
                                            SurfaceMeshWithFeatures& feature_surface)
{
    auto& m = feature_surface.surface_mesh;

    mesh::write_vtk(m, "debug_feature_surface.vtk");

    // Count # of interfaces
    mesh::ComponentIndex n_interfaces = 0;
    for (auto e : m.elements())
    {
        if (e.component_id() >= n_interfaces)
        {
            n_interfaces = e.component_id() + 1;
        }
    }

    // Create a map from each interface to the triangles that compose them
    std::vector<std::vector<mesh::LocalIndex>> interface_id_to_triangles_id_map(n_interfaces);
    for (auto e : m.elements())
    {
        interface_id_to_triangles_id_map[e.component_id()].push_back(e.local_index());
    }

    // Fill this structure with the PLC information
    tetgen::tetgenio in;
    in.firstnumber = 0; // Indicate indices should be zero-based

    // Populate the facets of the PLC
    int n_tetgen_facets = 0;
    std::vector<bool> is_planar_by_interface(interface_id_to_triangles_id_map.size(), false);
    for (auto interface_id : utility::indices(interface_id_to_triangles_id_map))
    {
        bool interface_is_planar = true;
        for (auto& global_facet_id :
             feature_surface.region_boundary_id_to_facet_ids_map[interface_id])
        {
            interface_is_planar &= is_facet_planar[global_facet_id];
        }

        if (interface_is_planar)
        {
            n_tetgen_facets += 1;
            is_planar_by_interface[interface_id] = true;
        }
        else
        {
            n_tetgen_facets += interface_id_to_triangles_id_map[interface_id].size();
        }
    }
    in.numberoffacets  = n_tetgen_facets;
    in.facetlist       = new tetgen::tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // Find the unmatched edges for each planar interface
    int curr_tetgen_facet_id     = 0;
    tetgen::tetgenio::facet* f   = nullptr;
    tetgen::tetgenio::polygon* p = nullptr;
    std::vector<mesh::LocalIndex> tetgen_to_mesh_point_index_map;
    std::unordered_map<mesh::LocalIndex, mesh::LocalIndex> mesh_to_tetgen_point_index_map;
    for (auto interface_id : utility::indices(interface_id_to_triangles_id_map))
    {
        bool interface_is_planar = true;
        for (auto& global_facet_id :
             feature_surface.region_boundary_id_to_facet_ids_map[interface_id])
        {
            interface_is_planar &= is_facet_planar[global_facet_id];
        }

        if (interface_is_planar)
        {
            std::unordered_map<std::size_t, std::pair<mesh::LocalIndex, mesh::LocalIndex>>
              unmatched_edges;
            for (auto& triangle_id : interface_id_to_triangles_id_map[interface_id])
            {
                auto e    = m.elements()[triangle_id];
                auto conn = e.connectivity();
                for (auto segment : mesh::polygon_segments(e.connectivity()))
                {
                    // check segment points to conn
                    auto hash = utility::hash_sorted(conn[segment.first], conn[segment.second]);
                    auto find_iter = unmatched_edges.find(hash);
                    if (find_iter == unmatched_edges.end())
                    {
                        unmatched_edges[hash] = { conn[segment.first], conn[segment.second] };
                    }
                    else
                    {
                        unmatched_edges.erase(find_iter);
                    }
                }
            }

            // Create structures to accelerate construction of the boundary loops of the
            // interface
            std::vector<std::array<mesh::LocalIndex, 2>> unmatched_edges_vector;
            for (auto& key_value : unmatched_edges)
            {
                unmatched_edges_vector.push_back(
                  { key_value.second.first, key_value.second.second });
            }


            // Key: node index
            // Value: { connected neighbor node index along unmatched edge for which this node
            //             index is first in the edge (edge going out from this node),
            //          connected neighbor node index along unmatched edge for which this node
            //             index is second in the edge (edge going into this node) }
            std::unordered_map<mesh::LocalIndex, std::array<mesh::LocalIndex, 2>>
              boundary_node_adjacency_list;
            for (auto i : utility::indices(unmatched_edges_vector))
            {
                assert(unmatched_edges_vector[i][0] != unmatched_edges_vector[i][1]);

                add_directed_edge(boundary_node_adjacency_list,
                                  unmatched_edges_vector[i][0],
                                  unmatched_edges_vector[i][1]);
            }

            for (auto& key_value : boundary_node_adjacency_list)
            {
                if (mesh_to_tetgen_point_index_map.count(key_value.first) == 0)
                {
                    auto new_index = tetgen_to_mesh_point_index_map.size();
                    tetgen_to_mesh_point_index_map.push_back(key_value.first);
                    mesh_to_tetgen_point_index_map[key_value.first] = new_index;
                }
            }

            // Build loops
            // Note: if this is a performance bottleneck using associative containers, renumber
            // the boundary nodes so that vectors can be used (keeping a map to the original
            // indices).
            std::unordered_set<mesh::LocalIndex> visited_nodes;
            std::vector<std::vector<mesh::LocalIndex>> boundary_loops;
            while (visited_nodes.size() != boundary_node_adjacency_list.size())
            {
                auto loop_start_iter =
                  std::find_if(boundary_node_adjacency_list.begin(),
                               boundary_node_adjacency_list.end(),
                               [&visited_nodes](const auto& key_value)
                               { return visited_nodes.count(key_value.first) == 0; });
                auto& start_node_index = loop_start_iter->first;
                // auto& end_node_index   = loop_start_iter->second[1];

                // Start a new loop
                boundary_loops.push_back({});

                mesh::LocalIndex curr_node_index = start_node_index;
                mesh::LocalIndex next_node_index = loop_start_iter->second[0];
                do
                {
                    boundary_loops.back().push_back(curr_node_index);
                    visited_nodes.insert(curr_node_index);
                    curr_node_index = next_node_index;
                    next_node_index = boundary_node_adjacency_list[curr_node_index][0];
                } while (curr_node_index != start_node_index);
            }

            // With the loop determine which direction each goes (true = outward, false =
            // inward)
            std::vector<bool> boundary_loops_outward_indicator(boundary_loops.size(), false);
            auto interface_normal = mesh::calc_polygon_unit_normal(
              m.elements()[interface_id_to_triangles_id_map[interface_id][0]]);
            int n_inverted_loops = 0;
            for (auto loop_index : utility::indices(boundary_loops))
            {
                std::vector<simple_math::Point> loop_points;
                loop_points.reserve(boundary_loops.size());
                for (auto& n : boundary_loops[loop_index])
                {
                    loop_points.push_back(m.nodes()[n].point());
                }
                auto loop_normal = mesh::calc_polygon_unit_normal(loop_points);
                if (interface_normal.dot(loop_normal) > 0)
                {
                    boundary_loops_outward_indicator[loop_index] = true;
                }
                else
                {
                    n_inverted_loops++;
                }
            }

            using namespace simple_math::point_operators;

            // Populate tetgen PLC
            in.facetmarkerlist[curr_tetgen_facet_id] = interface_id;
            f                                        = &in.facetlist[curr_tetgen_facet_id];
            f->numberofpolygons                      = (int)boundary_loops.size();
            f->polygonlist      = new tetgen::tetgenio::polygon[f->numberofpolygons];
            f->numberofholes    = n_inverted_loops;
            f->holelist         = new double[3 * n_inverted_loops];
            int curr_hole_index = 0;
            for (auto loop_index : utility::indices(boundary_loops))
            {
                // If the normal vector of the polygon formed by the loop points inward,
                // then it needs to be marked as a hole
                if (!boundary_loops_outward_indicator[loop_index])
                {
                    simple_math::Point hole_coord = { 0, 0, 0 };
                    for (auto node_index : boundary_loops[loop_index])
                    {
                        hole_coord = hole_coord + m.nodes()[node_index].point();
                    }
                    hole_coord = hole_coord / (double)boundary_loops[loop_index].size();
                    f->holelist[3 * curr_hole_index]     = hole_coord[0];
                    f->holelist[3 * curr_hole_index + 1] = hole_coord[1];
                    f->holelist[3 * curr_hole_index + 2] = hole_coord[2];
                    curr_hole_index++;
                }

                auto& loop          = boundary_loops[loop_index];
                p                   = &f->polygonlist[loop_index];
                p->numberofvertices = loop.size();
                p->vertexlist       = new int[p->numberofvertices];
                for (auto i : utility::range(loop.size()))
                {
                    p->vertexlist[i] = loop[i];
                }
            }
            curr_tetgen_facet_id++;
        } // end if the interface is planar
        else
        { // a non-planar interface
            auto triangle_ids = interface_id_to_triangles_id_map[interface_id];

            // Populate tetgen PLC
            // Since the interface is not planar, each triangle should be inserted as a facet
            for (auto& triangle_id : triangle_ids)
            {
                in.facetmarkerlist[curr_tetgen_facet_id] = interface_id;
                f                                        = &in.facetlist[curr_tetgen_facet_id];
                f->numberofpolygons                      = 1;
                f->polygonlist   = new tetgen::tetgenio::polygon[f->numberofpolygons];
                f->numberofholes = 0;
                f->holelist      = nullptr;

                p                   = &f->polygonlist[0];
                p->numberofvertices = 3;
                p->vertexlist       = new int[p->numberofvertices];
                auto triangle       = m.elements()[triangle_id];
                for (auto i : utility::range(3))
                {
                    auto node_id     = triangle.connectivity()[i];
                    p->vertexlist[i] = node_id;

                    if (mesh_to_tetgen_point_index_map.count(node_id) == 0)
                    {
                        auto new_index = tetgen_to_mesh_point_index_map.size();
                        tetgen_to_mesh_point_index_map.push_back(node_id);
                        mesh_to_tetgen_point_index_map[node_id] = new_index;
                    }
                }

                curr_tetgen_facet_id++;
            }
        }
    }


    // Remap the node indices to new numbering for tetgen that is 0-based and contiguous
    for (auto facet_id : utility::range(in.numberoffacets))
    {
        auto& f = in.facetlist[facet_id];
        for (auto polygon_id : utility::range(f.numberofpolygons))
        {
            auto& p = f.polygonlist[polygon_id];
            for (auto vertex_id : utility::range(p.numberofvertices))
            {
                p.vertexlist[vertex_id] = mesh_to_tetgen_point_index_map[p.vertexlist[vertex_id]];
            }
        }
    }

    // Populate the points of the PLC, now that we know which points are needed for the facets
    in.numberofpoints = static_cast<int>(tetgen_to_mesh_point_index_map.size());
    in.pointlist      = new REAL[in.numberofpoints * 3];
    for (auto tetgen_id : utility::indices(tetgen_to_mesh_point_index_map))
    {
        auto n = m.nodes()[tetgen_to_mesh_point_index_map[tetgen_id]];
        for (auto i : utility::range(3))
        {
            in.pointlist[3 * tetgen_id + i] = n.point(i);
        }
    }

    /*
    // Populate the points of the PLC, now that we know which points are needed for the facets
    in.numberofpoints = static_cast<int>(m.nodes().size());
    in.pointlist      = new REAL[in.numberofpoints * 3];
    for (auto n : m.nodes())
    {
        for (auto i : utility::range(3))
        {
            in.pointlist[3 * n.local_index() + i] = n.point(i);
        }
    }
    */
    // mesh::write_vtk(tmp_holes, "debug_holes.vtk");

    // Create tetrahedral mesh structure
    tetgen::tetgenbehavior flags;
    // Switches
    flags.plc    = 1;
    flags.object = tetgen::tetgenbehavior::POLY;
    // flags.nobisect    = 0;
    flags.quality = 1;
    flags.quiet   = 0;
    // flags.maxvolume   = 0.05;
    flags.mindihedral = 1.0;
    // flags.minratio    = 1.3;

    std::unique_ptr<tetgen::tetgenio> out_ptr = std::make_unique<tetgen::tetgenio>();
    tetgen::tetgenio& out                     = *out_ptr;

    {
        auto m = mesh::get_mesh_from_tetgen_input(in);
        mesh::write_vtk(m, "tetgen_input.vtk");
    }

    // The external call
    tetgen::tetrahedralize(&flags, &in, &out);

    /*
    flags.refine = 1;
    flags.quality     = 1;
    std::unique_ptr<tetgen::tetgenio> out_ptr2 = std::make_unique<tetgen::tetgenio>();
    tetgen::tetgenio& out2                     = *out_ptr2;
    tetgen::tetrahedralize(&flags, &out, &out2);
    */

    return mesh::get_mesh_from_tetgen_output(out);
}

#pragma endregion

#pragma region approach_3

PLC feature_surface_to_plc(const SurfaceMeshWithFeatures& feature_surface,
                           const std::vector<bool>& is_facet_planar)
{
    auto& m = feature_surface.surface_mesh;

    // Count # of interfaces
    mesh::ComponentIndex n_interfaces = 0;
    for (auto e : m.elements())
    {
        if (e.component_id() >= n_interfaces)
        {
            n_interfaces = e.component_id() + 1;
        }
    }

    // Create a map from each interface to the triangles that compose them
    std::vector<std::vector<mesh::LocalIndex>> interface_id_to_triangles_id_map(n_interfaces);
    for (auto e : m.elements())
    {
        interface_id_to_triangles_id_map[e.component_id()].push_back(e.local_index());
    }

    PLC plc;
    for (auto interface_id : utility::indices(interface_id_to_triangles_id_map))
    {
        bool interface_is_planar = true;
        for (auto& global_facet_id :
             feature_surface.region_boundary_id_to_facet_ids_map[interface_id])
        {
            interface_is_planar &= is_facet_planar[global_facet_id];
        }

        if (interface_is_planar)
        {
            std::unordered_map<std::size_t, std::pair<mesh::LocalIndex, mesh::LocalIndex>>
              unmatched_edges;
            for (auto& triangle_id : interface_id_to_triangles_id_map[interface_id])
            {
                auto e    = m.elements()[triangle_id];
                auto conn = e.connectivity();
                for (auto segment : mesh::polygon_segments(e.connectivity()))
                {
                    // check segment points to conn
                    auto hash = utility::hash_sorted(conn[segment.first], conn[segment.second]);
                    auto find_iter = unmatched_edges.find(hash);
                    if (find_iter == unmatched_edges.end())
                    {
                        unmatched_edges[hash] = { conn[segment.first], conn[segment.second] };
                    }
                    else
                    {
                        unmatched_edges.erase(find_iter);
                    }
                }
            }

            // Create structures to accelerate construction of the boundary loops of the
            // interface
            std::vector<std::array<mesh::LocalIndex, 2>> unmatched_edges_vector;
            for (auto& key_value : unmatched_edges)
            {
                unmatched_edges_vector.push_back(
                  { key_value.second.first, key_value.second.second });
            }

            // Key: node index
            // Value: { connected neighbor node index along unmatched edge for which this node
            //             index is first in the edge (edge going out from this node),
            //          connected neighbor node index along unmatched edge for which this node
            //             index is second in the edge (edge going into this node) }
            std::unordered_map<mesh::LocalIndex, std::array<mesh::LocalIndex, 2>>
              boundary_node_adjacency_list;
            for (auto i : utility::indices(unmatched_edges_vector))
            {
                assert(unmatched_edges_vector[i][0] != unmatched_edges_vector[i][1]);

                add_directed_edge(boundary_node_adjacency_list,
                                  unmatched_edges_vector[i][0],
                                  unmatched_edges_vector[i][1]);
            }

            // Populate the facet
            PLCFacet facet;

            // Build loops
            // Note: if this is a performance bottleneck using associative containers, renumber
            // the boundary nodes so that vectors can be used (keeping a map to the original
            // indices).
            std::unordered_set<mesh::LocalIndex> visited_nodes;
            while (visited_nodes.size() != boundary_node_adjacency_list.size())
            {
                auto loop_start_iter =
                  std::find_if(boundary_node_adjacency_list.begin(),
                               boundary_node_adjacency_list.end(),
                               [&visited_nodes](const auto& key_value)
                               { return visited_nodes.count(key_value.first) == 0; });
                auto& start_node_index = loop_start_iter->first;
                // auto& end_node_index   = loop_start_iter->second[1];

                // Start a new loop
                facet.boundary_loops.push_back({});

                mesh::LocalIndex curr_node_index = start_node_index;
                mesh::LocalIndex next_node_index = loop_start_iter->second[0];
                do
                {
                    facet.boundary_loops.back().push_back(curr_node_index);
                    visited_nodes.insert(curr_node_index);
                    curr_node_index = next_node_index;
                    next_node_index = boundary_node_adjacency_list[curr_node_index][0];
                } while (curr_node_index != start_node_index);
            }

            // With the loop determine which direction each goes (true = outward, false =
            // inward)
            facet.is_boundary_loop_a_hole.resize(facet.boundary_loops.size());
            auto interface_normal = mesh::calc_polygon_unit_normal(
              m.elements()[interface_id_to_triangles_id_map[interface_id][0]]);
            for (auto loop_index : utility::indices(facet.boundary_loops))
            {
                std::vector<simple_math::Point> loop_points;
                loop_points.reserve(facet.boundary_loops[loop_index].size());
                for (auto& n : facet.boundary_loops[loop_index])
                {
                    loop_points.push_back(m.nodes()[n].point());
                }

                auto loop_normal = mesh::calc_polygon_unit_normal(loop_points);
                if (interface_normal.dot(loop_normal) > 0)
                {
                    facet.is_boundary_loop_a_hole[loop_index] = true;
                }
                else
                {
                    facet.is_boundary_loop_a_hole[loop_index] = false;
                }
            }

            plc.facets.push_back(facet);
            // End populating the facet
        }
    }

    return plc;
}

// Triangulate each facet of the PLC
void plc_to_polygonal_surface(const PLC& plc,
                              const SurfaceMeshWithFeatures& feature_surface,
                              const double max_triangle_area)
{
    // Find the connectivity between loops (which facets share them)
    // If two facets share an entire loop, then one must be a hole.

    // How do I deal with the shared segments (partial loops)?

    // Loop over the facets
    // 1) Create a triangle mesh for the facet
    // 2) Update adjacent facets to refine shared partial/full loops

    using namespace simple_math::point_operators;

    // To prove the concept, mesh each independently
    auto& m = feature_surface.surface_mesh;
    for (auto facet_id : utility::indices(plc.facets))
    {
        auto& facet = plc.facets[facet_id];

        // Define a normalized coordinate system that is 2D baed on the plane points of the
        // first loop, which should be the same plane as other loops in this facet.
        simple_math::Matrix3D transformation;
        {
            std::vector<simple_math::Point> loop_points;
            simple_math::Point mean_point = { 0, 0, 0 };
            loop_points.reserve(facet.boundary_loops[0].size());
            for (auto& n : facet.boundary_loops[0])
            {
                loop_points.push_back(m.nodes()[n].point());
                mean_point = mean_point + loop_points.back();
            }
            mean_point = mean_point / (double)loop_points.size();

            auto axis_1 =
              simple_math::point_as_vector(loop_points.front() - mean_point).normalized();
            auto axis_3 = mesh::calc_polygon_unit_normal(loop_points);
            auto axis_2 = axis_3.cross(axis_1);

            transformation << axis_1[0], axis_1[1], axis_1[2], axis_2[0], axis_2[1], axis_2[2],
              axis_3[0], axis_3[1], axis_3[2];
        }

        // Find the subset of points needed for this facet
        std::vector<mesh::LocalIndex> facet_to_mesh_point_index_map;
        std::unordered_map<mesh::LocalIndex, mesh::LocalIndex> mesh_to_facet_point_index_map;
        for (auto& loop : facet.boundary_loops)
        {
            for (auto& node_id : loop)
            {
                if (mesh_to_facet_point_index_map.count(node_id) == 0)
                {
                    auto new_index = facet_to_mesh_point_index_map.size();
                    facet_to_mesh_point_index_map.push_back(node_id);
                    mesh_to_facet_point_index_map[node_id] = new_index;
                }
            }
        }

        // Transform coordinates to this normalized coordinate system
        std::vector<mesh::Point> transformed_points;
        // Save the nodes' indices in the connectivity as markers so we can recover the nodes
        // after triangulation
        std::vector<int> markers;
        double normalized_z_offset = 0.0;
        mesh::LocalIndex count     = 0;
        for (auto& mesh_node_id : facet_to_mesh_point_index_map)
        {
            auto transformed_p  = transformation * m.nodes()[mesh_node_id].coord();
            normalized_z_offset = transformed_p(2);

            transformed_points.push_back({ transformed_p(0), transformed_p(1), 0.0 });
            markers.push_back(count++);
        }

        // Create the boundary segments for the triangulation
        std::vector<std::pair<mesh::LocalIndex, mesh::LocalIndex>> boundary_segments;
        for (auto& loop : facet.boundary_loops)
        {
            for (auto s : mesh::polygon_segments(
                   utility::Slice<const mesh::LocalIndex, mesh::ConnIndex>(loop.data(),
                                                                           loop.size())))
            {
                boundary_segments.push_back({ mesh_to_facet_point_index_map[s.first],
                                              mesh_to_facet_point_index_map[s.second] });
            }
        }

        // Call triangle
        triangulateio output;
        try
        {
            output = mesh::detail::triangle_polygon_kernel(transformed_points,
                                                           markers,
                                                           boundary_segments,
                                                           {},
                                                           {},
                                                           {},
                                                           {},
                                                           false,
                                                           max_triangle_area,
                                                           10.0);

            mesh::Mesh tmp;
            for (auto& p : transformed_points)
            {
                tmp.nodes().add(p);
            }
            for (auto& s : boundary_segments)
            {
                tmp.elements().add(mesh::ElementType::L1D2N, { s.first, s.second });
            }
            write_vtk(tmp, "tmp.vtk", false);
        }
        catch (...)
        {
            throw _UTIL_MSG_EXCEPT("Triangle failed to mesh the facet representing an interface.  "
                                   "TODO output a VTK of the failed input.");
        }

        {
            // So, I will exploit the fact that triangle labels the given points first
            // Add all the points, starting after the points I seeded
            transformation.transposeInPlace(); // Now this transforms normalized to global
            mesh::Mesh triangulation;
            for (auto i : utility::range(output.numberofpoints))
            {
                // Note: we need to transform the coords back to global system
                simple_math::Vector3D coord_2D = { output.pointlist[i * 2],
                                                   output.pointlist[i * 2 + 1],
                                                   normalized_z_offset };

                auto global_coord = transformation * coord_2D;
                if (!global_coord.array().isFinite().all())
                {
                    // See mesh::trianulate_polygon for reference of error info
                    throw _UTIL_MSG_EXCEPT("Transformation failed.");
                }
                triangulation.nodes().add(simple_math::to_point(global_coord));
            }
            // Now add the triangles
            std::vector<mesh::LocalIndex> conn(3, 0);
            for (int i = 0; i < output.numberoftriangles; i++)
            {
                for (auto j : utility::range(3))
                {
                    conn[j] = output.trianglelist[i * 3 + j];
                }
                triangulation.elements().add(mesh::ElementType::T2D3N, conn);
                triangulation.elements().back().component_id() = facet_id;
            }
            mesh::write_vtk(triangulation, "facet_" + std::to_string(facet_id) + ".vtk");
        }

        // frees memory used by triangle
        free(output.pointlist);
        free(output.pointattributelist);
        free(output.pointmarkerlist);
        free(output.trianglelist);
        free(output.triangleattributelist);
        free(output.neighborlist);
        free(output.segmentlist);
        free(output.segmentmarkerlist);
        free(output.edgelist);
        free(output.edgemarkerlist);
    }
}

void polygonal_surface_to_tetrahedralization()
{
    //
}


#pragma endregion

#endif
}