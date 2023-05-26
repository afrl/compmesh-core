#include "comp_mesh_core/FiberBundle.h"

#include "mesh/Mesh.h"
#include "mesh/algorithms/IO.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/StringUtils.h"

namespace composite_mesher
{
FiberBundle read_fib_file(const std::string& file_path)
{
    // Open the ctw file
    std::ifstream fin(file_path.c_str());
    if (!fin)
    {
        throw _UTIL_MSG_EXCEPT("Cannot find file " + file_path);
    }

    utility::skip_to_line_that_contains("*FIBERS", fin);

    std::string line;
    // Expects: "#Fiber	#Material	#Diameter"
    utility::get_next_line(fin, line);

    // Expects: <# of fibers> <material id> <diameter>
    utility::get_next_line(fin, line);
    auto tokens = utility::tokenize(line);

    FiberBundle bundle;
    bundle.fibers.resize(atoi(tokens[0].c_str()));
    bundle.material_id = atoi(tokens[1].c_str());
    bundle.diameter    = atof(tokens[2].c_str());

    for (auto i : utility::indices(bundle.fibers))
    {
        // Expects: <fiber index> <# of points>
        utility::get_next_line(fin, line);
        tokens = utility::tokenize(line);
        bundle.fibers[i].points.resize(atoi(tokens[1].c_str()));
    }

    // Expects: "Node	X	Y	Z" for 1st (blank lines after 1st section)
    utility::get_next_line(fin, line);

    for (auto i : utility::indices(bundle.fibers))
    {
        // Expects: <fiber index>
        utility::get_next_line(fin, line);

        for (auto j : utility::indices(bundle.fibers[i].points))
        {
            // Expects: <point index (1-based)> <x> <y> <z>
            utility::get_next_line(fin, line);
            tokens  = utility::tokenize(line);
            auto& p = bundle.fibers[i].points[j];
            for (auto k : utility::range(3))
            {
                p[k] = atof(tokens[k + 1].c_str());
            }
        }
    }

    return bundle;
}

mesh::Mesh to_mesh(const FiberBundle& bundle)
{
    mesh::Mesh bundle_mesh;
    for (auto i : utility::indices(bundle.fibers))
    {
        auto offset = bundle_mesh.nodes().size();
        for (auto j : utility::indices(bundle.fibers[i].points))
        {
            bundle_mesh.nodes().add(bundle.fibers[i].points[j], j);
        }
        for (auto j : utility::range(bundle.fibers[i].points.size() - 1))
        {
            bundle_mesh.elements().add(mesh::ElementType::L1D2N,
                                      { offset + mesh::LocalIndex(j), offset + mesh::LocalIndex(j) + 1 },
                                      j);
            bundle_mesh.elements().back().component_id() = i;
        }
    }
    return bundle_mesh;
}

PointMatrix extract_points_matrix(const FiberBundle& bundle)
{
    int n_points = 0;
    for (auto& fiber : bundle.fibers)
    {
        n_points += fiber.points.size();
    }
    Eigen::MatrixX3d points = Eigen::MatrixX3d::Zero(n_points, 3);
    std::vector<Index> fiber_ids(n_points);
    std::vector<Index> index_in_fiber(n_points);

    int curr_index = 0;
    for (auto i : utility::indices(bundle.fibers))
    {
        auto& fiber = bundle.fibers[i];
        for (auto j : utility::indices(fiber.points))
        {
            auto& p = fiber.points[j];
            for (auto k : utility::range(3))
            {
                points(curr_index, k) = p[k];
            }
            fiber_ids[curr_index]      = i;
            index_in_fiber[curr_index] = j;
            curr_index++;
        }
    }

    return { points, fiber_ids, index_in_fiber };
}

FiberDistanceInfo get_distance_info(Eigen::Map<const simple_math::DynamicDoubleMatrix> points_ref)
{
    simple_math::DynamicVector<double> distance_along_line(points_ref.rows());
    distance_along_line.setZero();
    for (auto i : utility::range(points_ref.rows() - 1))
    {
        const double d                = (points_ref.row(i + 1) - points_ref.row(i)).norm();
        distance_along_line(i + 1, 0) = distance_along_line(i, 0) + d;
    }
    double total_length = distance_along_line(points_ref.rows() - 1);

    return { total_length, distance_along_line, distance_along_line / total_length };
}

double compute_distance_along_fiber(const Index& intersected_segment,
                                    const double& segment_intersection_fraction,
                                    const FiberDistanceInfo& distance_info)
{
    const Index& segment_p_start = intersected_segment;
    const Index& segment_p_end   = intersected_segment + 1;

    // Compute the actual distance along this fiber for the last cross-section
    const auto& p_start_distance = distance_info.distance_along_line[segment_p_start];
    const auto& p_end_distance   = distance_info.distance_along_line[segment_p_end];
    const double distance =
      segment_intersection_fraction * (p_end_distance - p_start_distance) + p_start_distance;
    return distance;
}

std::vector<double> compute_distance_along_fibers(
  const std::vector<Index>& intersected_segments,
  const std::vector<double>& segment_intersection_fractions,
  const std::vector<FiberDistanceInfo>& distance_infos)
{
    std::vector<double> distances(distance_infos.size());
    for (auto fiber_id : utility::indices(distance_infos))
    {
        distances[fiber_id] = compute_distance_along_fiber(intersected_segments[fiber_id],
                                                           segment_intersection_fractions[fiber_id],
                                                           distance_infos[fiber_id]);
    }
    return distances;
}

std::vector<double> compute_fractions_along_fibers(
  const std::vector<Index>& intersected_segments,
  const std::vector<double>& segment_intersection_fractions,
  const std::vector<FiberDistanceInfo>& distance_infos)
{
    std::vector<double> normalized_distances(distance_infos.size());
    for (auto fiber_id : utility::indices(distance_infos))
    {
        normalized_distances[fiber_id] =
          compute_distance_along_fiber(intersected_segments[fiber_id],
                                       segment_intersection_fractions[fiber_id],
                                       distance_infos[fiber_id]) /
          distance_infos[fiber_id].total_length;
    }
    return normalized_distances;
}
}