#include "comp_mesh_core/TowSurfaceGeneration.h"

#include "comp_mesh_core/GlobalOptions.h"
#include "mesh/algorithms/IO.h"
#include "simple_math/ConvexHull.h"
#include "simple_math/PlaneGeometry.h"
#include "simple_math/SplineLine.h"
#include "simple_math/Utility.h"

#include <cmath>
#include <vector>

namespace composite_mesher
{

VTMSTowSurface generate_surface(const FiberBundle& bundle,
                                const int nominal_n_cross_sections,
                                const int n_points_around_each_cross_section,
                                const double max_angle_allowed_between_fibers_and_plane,
                                const double min_angle_allowed_between_cross_sections,
                                const int n_cutbacks_allowed)
{
    auto cross_section_planes_and_intersections =
      identify_cross_section_planes_and_intersections(bundle,
                                                      nominal_n_cross_sections,
                                                      max_angle_allowed_between_fibers_and_plane,
                                                      min_angle_allowed_between_cross_sections,
                                                      n_cutbacks_allowed);

    auto discrete_cross_sections = discretize_cross_sections(cross_section_planes_and_intersections,
                                                             n_points_around_each_cross_section);


    if (GlobalOptions::get().write_intermediate_meshes)
    {
        mesh::Mesh intersections;
        std::vector<uint32_t> is_convex_hull;
        std::vector<uint32_t> intersected_segments;
        std::vector<double> intersection_lengths;
        for (auto i : utility::indices(discrete_cross_sections))
        {
            auto& cs_plane = cross_section_planes_and_intersections[i];
            auto& cs       = discrete_cross_sections[i];
            for (auto j : utility::indices(cs.fiber_intersections))
            {
                auto& inter = cs.fiber_intersections[j];
                intersections.nodes().add(simple_math::to_point(inter), i);
                is_convex_hull.push_back(cs.is_fiber_convex_hull[j]);
                intersected_segments.push_back(cs.fiber_intersected_segments[j]);
                intersection_lengths.push_back(cs_plane.fiber_intersection_length[j]);
            }
        }
        mesh::MeshAssociatedData data;
        data.add_node_uint32("is_convex_hull", is_convex_hull);
        data.add_node_uint32("intersected_segments", intersected_segments);
        data.add_node_double("intersection_lengths", intersection_lengths);
        mesh::write_vtk(intersections, "tow_intersects.vtk", data);


        mesh::Mesh fit_cross_section;
        for (auto i : utility::indices(discrete_cross_sections))
        {
            auto& cs = discrete_cross_sections[i];
            for (auto p : cs.discrete_cross_section_points)
            {
                fit_cross_section.nodes().add(simple_math::to_point(p), i);
            }
        }
        mesh::write_vtk(fit_cross_section, "tow_cross_sections.vtk");
    }

    VTMSTowSurface tow_surface;
    tow_surface.surface_normal_direction = VTMSTowSurface::NormalDirection::OUTWARD;
    tow_surface.type                     = VTMSTowSurface::MeshType::POLY;
    for (auto& cs : discrete_cross_sections)
    {
        for (auto& p : cs.discrete_cross_section_points)
        {
            tow_surface.mesh.nodes().add(simple_math::to_point(p));
        }
        tow_surface.medials.push_back(simple_math::to_point(cs.point));
    }
    for (auto i : utility::range(std::size_t(1), discrete_cross_sections.size()))
    {
        auto& cs_1 = discrete_cross_sections[i - 1].discrete_cross_section_points;
        auto& cs_2 = discrete_cross_sections[i].discrete_cross_section_points;

        double min_distance     = std::numeric_limits<double>::max();
        int min_distance_offset = 0;
        for (auto j : utility::range(n_points_around_each_cross_section))
        {
            double distance_sum = 0.;
            for (auto k : utility::range(n_points_around_each_cross_section))
            {
                distance_sum +=
                  (cs_1[k] - cs_2[(k + j) % n_points_around_each_cross_section]).norm();
            }
            if (distance_sum < min_distance)
            {
                min_distance        = distance_sum;
                min_distance_offset = j;
            }
        }

        for (auto k : utility::range(n_points_around_each_cross_section))
        {
            // cross-section i-1    cs i
            //             p1        p3
            //             |         |
            //             a1        a2
            //             |         |
            //             p2        p4
            //
            // There are two possibilities:
            //
            // 1) cross-section i-1  cs i
            //                p1 --- p3
            //                 | \   |
            //                 |  \  |
            //                 |   \ |
            //                p2 --- p4
            //
            // 2) cross-section i-1  cs i
            //                p1 --- p3
            //                 |   / |
            //                 |  /  |
            //                 | /   |
            //                p2 --- p4
            //
            // We will pick which one based on the which option maximizes convexity.  To
            // determine which option results in convexity, p1->p2 crossed with p1->p4  should
            // be in the opposite direction of a1->a2 (where a1 is the average of p1 and p2 and
            // a2 is the average of p3 and p4).
            //                p1 --- p3
            //                 |     |
            //                a1 --> a2
            //                 |     |
            //                p2 --- p4

            Index p1_i = n_points_around_each_cross_section * (i - 1) + k;
            Index p2_i = n_points_around_each_cross_section * (i - 1) +
                         ((k + 1) % n_points_around_each_cross_section);
            Index p3_i = n_points_around_each_cross_section * i +
                         ((k + min_distance_offset) % n_points_around_each_cross_section);
            Index p4_i = n_points_around_each_cross_section * i +
                         ((k + min_distance_offset + 1) % n_points_around_each_cross_section);

            auto p1 = tow_surface.mesh.nodes()[p1_i].coord();
            auto p2 = tow_surface.mesh.nodes()[p2_i].coord();
            auto p3 = tow_surface.mesh.nodes()[p3_i].coord();
            auto p4 = tow_surface.mesh.nodes()[p4_i].coord();

            if ((p2 - p1).cross(p4 - p1).dot((p4 + p3) / 2. - (p2 + p1) / 2.) < 0)
            {
                tow_surface.mesh.elements().add(mesh::ElementType::T2D3N, { p1_i, p2_i, p4_i });
                tow_surface.mesh.elements().add(mesh::ElementType::T2D3N, { p1_i, p4_i, p3_i });
            }
            else
            {
                tow_surface.mesh.elements().add(mesh::ElementType::T2D3N, { p1_i, p2_i, p3_i });
                tow_surface.mesh.elements().add(mesh::ElementType::T2D3N, { p2_i, p4_i, p3_i });
            }
        }
    }
    tow_surface.surface_element_types.resize(tow_surface.mesh.elements().size(),
                                             VTMSTowSurface::SurfaceElementType::INTERIOR);

    // Add first end cap, note that order of cross-section points are oriented in the direction of
    // the tow path.
    std::vector<Index> cap_connectivity(
      discrete_cross_sections[0].discrete_cross_section_points.size());
    std::iota(cap_connectivity.begin(), cap_connectivity.end(), 0);
    // Reverse to make point outwards
    std::reverse(cap_connectivity.begin(), cap_connectivity.end());
    tow_surface.mesh.elements().add(mesh::ElementType::POLYGON, cap_connectivity);
    // And second end cap
    cap_connectivity.resize(discrete_cross_sections.back().discrete_cross_section_points.size());
    std::iota(cap_connectivity.begin(),
              cap_connectivity.end(),
              tow_surface.mesh.nodes().size() - cap_connectivity.size());
    tow_surface.mesh.elements().add(mesh::ElementType::POLYGON, cap_connectivity);
    tow_surface.surface_element_types.resize(tow_surface.mesh.elements().size(),
                                             VTMSTowSurface::SurfaceElementType::BOUNDARY);

    triangulate_tow_surface(tow_surface, false);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        mesh::write_vtk(tow_surface.mesh, "tow.vtk");
    }

    return tow_surface;
}

CrossSectionPlane identify_cross_sectional_plane(
  const FiberBundle& bundle,
  const std::vector<Eigen::Map<const simple_math::DynamicDoubleMatrix>>& fiber_point_matrices,
  const std::vector<FiberDistanceInfo>& distance_infos,
  const double& delta_s,
  const std::vector<Index>& prev_fiber_intersected_segments,
  const std::vector<double>& prev_fiber_segment_intersection_fractions)
{
    // For a relative increment (defined by delta_s) along each fiber, find the point and
    // orientation.  The average point and orientation of all fibers will define the plane for
    // the next cross-section.
    simple_math::Vector3D average_point                  = simple_math::Vector3D::Zero();
    simple_math::Vector3D average_tangent                = simple_math::Vector3D::Zero();
    std::vector<Index> extrap_fiber_intersected_segments = prev_fiber_intersected_segments;
    std::vector<Index> fiber_segments_trust_region_start = prev_fiber_intersected_segments;
    std::vector<Index> fiber_segments_trust_region_end   = prev_fiber_intersected_segments;
    std::vector<simple_math::Vector3D> fiber_tangents(bundle.fibers.size());
    for (auto fiber_id : utility::indices(bundle.fibers))
    {
        // Compute the actual distance along this fiber for the last cross-section
        const double prev_d =
          compute_distance_along_fiber(prev_fiber_intersected_segments[fiber_id],
                                       prev_fiber_segment_intersection_fractions[fiber_id],
                                       distance_infos[fiber_id]);

        // Compute the estimated distance along this fiber for the next cross-section
        const double delta_d              = delta_s * distance_infos[fiber_id].total_length;
        const double estimated_next_d     = prev_d + delta_d;
        const double trust_region_start_d = prev_d + delta_d / 2.0;
        const double trust_region_end_d   = prev_d + 2.0 * delta_d;

        // Find the segment that straddles the estimated distance for the next cross-section
        Index next_segment_p_end = prev_fiber_intersected_segments[fiber_id] + 1;
        double next_p_end_d      = distance_infos[fiber_id].distance_along_line[next_segment_p_end];
        if (distance_infos[fiber_id].total_length > estimated_next_d)
        {
            // Iterate along the fiber until the segment that straddles "trust_region_start_d" is
            // found, which will be the start of the trust region.
            while (next_p_end_d <= estimated_next_d / 2.0 &&
                   extrap_fiber_intersected_segments[fiber_id] !=
                     bundle.fibers[fiber_id].points.size() - 2)
            {
                extrap_fiber_intersected_segments[fiber_id]++;
                fiber_segments_trust_region_start[fiber_id]++;
                fiber_segments_trust_region_end[fiber_id]++;
                next_segment_p_end = extrap_fiber_intersected_segments[fiber_id] + 1;
                next_p_end_d = distance_infos[fiber_id].distance_along_line[next_segment_p_end];
            }
            // Iterate along the fiber until the segment that straddles "estimated_next_d" is found,
            // which will be the best guess for the search in
            // find_intersections_with_cross_sectional_plane()
            while (next_p_end_d <= estimated_next_d &&
                   extrap_fiber_intersected_segments[fiber_id] !=
                     bundle.fibers[fiber_id].points.size() - 2)
            {
                extrap_fiber_intersected_segments[fiber_id]++;
                fiber_segments_trust_region_end[fiber_id]++;
                next_segment_p_end = extrap_fiber_intersected_segments[fiber_id] + 1;
                next_p_end_d = distance_infos[fiber_id].distance_along_line[next_segment_p_end];
            }
            // Continue iterating the trust region end along the fiber until the segment that
            // straddles "trust_region_end_d" is found.
            while (distance_infos[fiber_id]
                       .distance_along_line[fiber_segments_trust_region_end[fiber_id] + 1] <=
                     trust_region_end_d &&
                   fiber_segments_trust_region_end[fiber_id] !=
                     bundle.fibers[fiber_id].points.size() - 2)
            {
                fiber_segments_trust_region_end[fiber_id]++;
            }
            if (fiber_segments_trust_region_start[fiber_id] ==
                  prev_fiber_intersected_segments[fiber_id] &&
                fiber_segments_trust_region_end[fiber_id] !=
                  prev_fiber_intersected_segments[fiber_id])
            {
                fiber_segments_trust_region_start[fiber_id]++;
            }
        }
        else
        {
            extrap_fiber_intersected_segments[fiber_id] = bundle.fibers[fiber_id].points.size() - 2;
            fiber_segments_trust_region_end[fiber_id] = extrap_fiber_intersected_segments[fiber_id];
            next_segment_p_end = extrap_fiber_intersected_segments[fiber_id] + 1;
            next_p_end_d       = distance_infos[fiber_id].distance_along_line[next_segment_p_end];
        }
        const Index next_segment_p_start = next_segment_p_end - 1;
        const double next_p_start_d =
          distance_infos[fiber_id].distance_along_line[next_segment_p_start];

        // Compute the point for the current segment and add it to the average point
        const double next_r = (estimated_next_d - next_p_start_d) / (next_p_end_d - next_p_start_d);
        const simple_math::Vector3D next_p =
          next_r * (fiber_point_matrices[fiber_id].row(next_segment_p_end) -
                    fiber_point_matrices[fiber_id].row(next_segment_p_start)) +
          fiber_point_matrices[fiber_id].row(next_segment_p_start);
        average_point += next_p;

        // Compute the tangent vector for the current segment and add it to the average tangent
        fiber_tangents.push_back((fiber_point_matrices[fiber_id].row(next_segment_p_end) -
                                  fiber_point_matrices[fiber_id].row(next_segment_p_start))
                                   .normalized());
        average_tangent += fiber_tangents.back();
    }
    average_point /= bundle.fibers.size();
    average_tangent /= bundle.fibers.size();
    average_tangent.normalize();

    return { average_point,
             average_tangent,
             fiber_segments_trust_region_start,
             fiber_segments_trust_region_end,
             extrap_fiber_intersected_segments,
             fiber_tangents };
}

CrossSectionIntersections find_intersections_with_cross_sectional_plane(
  const FiberBundle& bundle,
  const std::vector<Eigen::Map<const simple_math::DynamicDoubleMatrix>>& fiber_point_matrices,
  const CrossSectionPlane& plane)
{
    // We have the average tangent and point to define the normal vector the for next
    // cross-sectional plane. Now, find the intersection of the plane with each fiber.
    // Note: curr_segment_by_fiber has been updated to be an estimate and should be close to
    // the intersection.
    std::vector<simple_math::Vector3D> fiber_intersections(bundle.fibers.size());
    std::vector<Index> fiber_intersected_segments = plane.extrapolated_fiber_intersected_segments;
    std::vector<double> fiber_segment_intersection_fractions(bundle.fibers.size(), 0.);
    for (auto fiber_id : utility::indices(bundle.fibers))
    {
        Index next_segment_p_start = fiber_intersected_segments[fiber_id];
        Index next_segment_p_end   = fiber_intersected_segments[fiber_id] + 1;

        auto next_segment_p_start_signed_dist = plane.normal.dot(
          fiber_point_matrices[fiber_id].row(next_segment_p_start).transpose() - plane.point);
        auto next_segment_p_end_signed_dist = plane.normal.dot(
          fiber_point_matrices[fiber_id].row(next_segment_p_end).transpose() - plane.point);

        auto next_segment_p_start_sign = simple_math::sign(next_segment_p_start_signed_dist);
        auto next_segment_p_end_sign   = simple_math::sign(next_segment_p_end_signed_dist);

        // Indicates if the current segment is intersected by the plane
        bool is_intersected = next_segment_p_start_sign == 0 || next_segment_p_end_sign == 0 ||
                              next_segment_p_end_sign != next_segment_p_start_sign;

        // Use this to monitor if we are oscillating around the same point
        Index last_segment_guess = fiber_intersected_segments[fiber_id];

        // Iterate until the specific segment intersected by the plane is found
        while (!is_intersected)
        {
            // If the end of the segment is closer to the plane than the start, then search further
            // along the fiber.
            if (abs(next_segment_p_end_signed_dist) < abs(next_segment_p_start_signed_dist))
            {
                // If we reached the end of the trust region, then just use that segment
                // to find the intersection with the plane.
                if (fiber_intersected_segments[fiber_id] ==
                    plane.fiber_segments_trust_region_end[fiber_id])
                {
                    break;
                }
                // Check if the search is oscillating between two segments that form a
                // local minima for the signed distance function.  For example, jumping
                // between segment 1 and 2 below:
                //
                // Plane
                //    |   o
                //    |  / segment 1
                //    | o
                //    |  \ segment 2
                //    |   o
                else if (fiber_intersected_segments[fiber_id] + 1 == last_segment_guess)
                {
                    break;
                }
                else
                {
                    last_segment_guess = fiber_intersected_segments[fiber_id];
                    fiber_intersected_segments[fiber_id]++;
                }
            }
            // If the end of the segment is farther from the plane than the start, then search
            // towards the start of the fiber.
            else
            {
                // If we reached the start of the trust region, then just use that segment
                // to find the intersection with the plane.
                if (fiber_intersected_segments[fiber_id] ==
                    plane.fiber_segments_trust_region_start[fiber_id])
                {
                    break;
                }
                // Check if the search is oscillating between two segments that form a
                // local minima for the signed distance function.  For example, jumping
                // between segment 1 and 2 below:
                //
                // Plane
                //    |   o
                //    |  / segment 1
                //    | o
                //    |  \ segment 2
                //    |   o
                else if (fiber_intersected_segments[fiber_id] - 1 == last_segment_guess)
                {
                    break;
                }
                else
                {
                    last_segment_guess = fiber_intersected_segments[fiber_id];
                    fiber_intersected_segments[fiber_id]--;
                }
            }

            next_segment_p_start = fiber_intersected_segments[fiber_id];
            next_segment_p_end   = fiber_intersected_segments[fiber_id] + 1;

            next_segment_p_start_signed_dist = plane.normal.dot(
              fiber_point_matrices[fiber_id].row(next_segment_p_start).transpose() - plane.point);
            next_segment_p_end_signed_dist = plane.normal.dot(
              fiber_point_matrices[fiber_id].row(next_segment_p_end).transpose() - plane.point);

            next_segment_p_start_sign = simple_math::sign(next_segment_p_start_signed_dist);
            next_segment_p_end_sign   = simple_math::sign(next_segment_p_end_signed_dist);

            is_intersected = next_segment_p_start_sign == 0 || next_segment_p_end_sign == 0 ||
                             next_segment_p_end_sign != next_segment_p_start_sign;
        }

        // Note: curr_segment_by_fiber was updated above
        // If a segment was found that is intersected by the plane, then it is straightforward to
        // find the point of intersection.  However, if there were oscillations or the fiber ends
        // near the plane (with no intersection), then we need to use a project of the closest point
        // of the segments onto the plane.  For this case, it should be noted that the "fraction"
        // along the segment does not correspond to the "intersection point".
        simple_math::Vector3D intersection_point;
        if (is_intersected)
        {
            // Now, update fiber_segment_intersection_fractions
            fiber_segment_intersection_fractions[fiber_id] =
              -next_segment_p_start_signed_dist /
              (next_segment_p_end_signed_dist - next_segment_p_start_signed_dist);

            // Find the intersection point
            intersection_point =
              fiber_segment_intersection_fractions[fiber_id] *
                (fiber_point_matrices[fiber_id].row(next_segment_p_end).transpose() -
                 fiber_point_matrices[fiber_id].row(next_segment_p_start).transpose()) +
              fiber_point_matrices[fiber_id].row(next_segment_p_start).transpose();
        }
        else
        {
            // The starting point of the segment was closer so use that projection
            simple_math::Vector3D point;
            if (abs(next_segment_p_start_signed_dist) < abs(next_segment_p_end_signed_dist))
            {
                fiber_segment_intersection_fractions[fiber_id] = 0.;
                point = fiber_point_matrices[fiber_id].row(next_segment_p_start);
            }
            else
            {
                fiber_segment_intersection_fractions[fiber_id] = 1.;
                point = fiber_point_matrices[fiber_id].row(next_segment_p_end);
            }
            intersection_point = point - plane.normal.dot(point - plane.point) * plane.normal;
        }

        bool is_inf = Eigen::isinf(intersection_point.array()).any();
        bool is_nan = Eigen::isnan(intersection_point.array()).any();

        if (is_inf || is_nan)
        {
            std::cout
              << "Error detected for intersection between a cross-sectional plane and fiber "
              << fiber_id << ":\n";
            std::cout << intersection_point << std::endl;
        }

        fiber_intersections[fiber_id] = intersection_point;
    }

    return { fiber_intersections,
             fiber_intersected_segments,
             fiber_segment_intersection_fractions };
}

std::vector<IntersectedCrossSection> identify_cross_section_planes_and_intersections(
  const FiberBundle& bundle,
  const int n_cross_sections,
  const double max_angle_deviation_for_cross_section_allowed,
  const double min_angle_between_cross_sections_allowed,
  const int n_cutbacks_allowed)
{
    std::vector<Eigen::Map<const simple_math::DynamicDoubleMatrix>> fiber_point_matrices;
    fiber_point_matrices.reserve(bundle.fibers.size());
    for (auto& fiber : bundle.fibers)
    {
        fiber_point_matrices.push_back(simple_math::points_as_matrix(fiber.points));
    }

    std::vector<FiberDistanceInfo> distance_infos(bundle.fibers.size());
    for (auto fiber_id : utility::indices(bundle.fibers))
    {
        distance_infos[fiber_id] = get_distance_info(fiber_point_matrices[fiber_id]);
    }

    // Stores index of the segment intersected by the current cross-section for each fiber
    std::vector<Index> prev_fiber_intersected_segments(bundle.fibers.size(), 0);
    // Stores fractional distance along the current segment that lies on the current cross-section
    // for each fiber
    std::vector<double> prev_fiber_segment_intersection_fractions(bundle.fibers.size(), 0.);

    std::vector<IntersectedCrossSection> cross_sections;
    cross_sections.reserve(n_cross_sections);
    bool reached_end_of_a_fiber = false;
    const double min_angle_allowed_radians =
      min_angle_between_cross_sections_allowed * EIGEN_PI / 180.0;
    const double min_cross_section_dot_product = cos(min_angle_allowed_radians);
    const double min_dev_dot_product =
      cos(max_angle_deviation_for_cross_section_allowed * EIGEN_PI / 180.0);
    while (!reached_end_of_a_fiber)
    {
        // For the first, cross-section, delta_s is 0, afterwards it is 1/(N-1), which will end at 1
        double delta_s = cross_sections.size() > 0 ? 1.0 / double(n_cross_sections - 1) : 0;

        CrossSectionPlane plane;
        bool angle_between_fibers_and_normal_criterion_triggered = false;
        bool angle_between_fibers_and_normal_criterion_satisfied = true;
        bool angle_between_cs_criterion_satisfied                = false;
        int n_bisections                                         = 0;
        do
        {
            plane = identify_cross_sectional_plane(bundle,
                                                   fiber_point_matrices,
                                                   distance_infos,
                                                   delta_s,
                                                   prev_fiber_intersected_segments,
                                                   prev_fiber_segment_intersection_fractions);

            if (cross_sections.size() > 0)
            {
                /*
                // If the fibers on the cross-section are highly unaligned, then increase increment
                // to get past the spurious point.
                angle_between_fibers_and_normal_criterion_satisfied = true;
                for (auto& tangent : plane.fiber_tangents)
                {
                    if (tangent.dot(plane.normal) < min_dev_dot_product)
                    {
                        angle_between_fibers_and_normal_criterion_satisfied = false;
                        angle_between_fibers_and_normal_criterion_triggered = true;
                        delta_s *= 2.0;
                        break;
                    }
                }
                */
                // If normal changes drastically then back off on increment.
                // Note: the above constraint overrides this one.
                if (!angle_between_fibers_and_normal_criterion_triggered)
                {
                    if (plane.normal.dot(cross_sections.back().normal) <
                        min_cross_section_dot_product)
                    {
                        delta_s /= 2.0;
                        n_bisections++;
                    }
                    else
                    {
                        angle_between_cs_criterion_satisfied = true;
                    }
                }
            }
        } while (!angle_between_fibers_and_normal_criterion_satisfied ||
                 (!angle_between_cs_criterion_satisfied && n_bisections < n_cutbacks_allowed));

        auto intersected_cs =
          find_intersections_with_cross_sectional_plane(bundle, fiber_point_matrices, plane);

        auto distances_along_fibers =
          compute_distance_along_fibers(intersected_cs.fiber_intersected_segments,
                                        intersected_cs.fiber_segment_intersection_fractions,
                                        distance_infos);

        cross_sections.push_back({ plane.point,
                                   plane.normal,
                                   intersected_cs.fiber_intersections,
                                   intersected_cs.fiber_intersected_segments,
                                   intersected_cs.fiber_segment_intersection_fractions,
                                   distances_along_fibers });

        prev_fiber_intersected_segments = intersected_cs.fiber_intersected_segments;
        prev_fiber_segment_intersection_fractions =
          intersected_cs.fiber_segment_intersection_fractions;

        auto normalized_distances =
          compute_fractions_along_fibers(prev_fiber_intersected_segments,
                                         prev_fiber_segment_intersection_fractions,
                                         distance_infos);
        if (*std::max_element(normalized_distances.begin(), normalized_distances.end()) >= 1.0)
        {
            reached_end_of_a_fiber = true;
        }
    }

    return cross_sections;
}

std::vector<DiscretizedCrossSection> discretize_cross_sections(
  const std::vector<IntersectedCrossSection>& cs_planes_and_intersections,
  const int n_points_around_each_cross_section)
{
    std::vector<DiscretizedCrossSection> discretized_cross_sections(
      cs_planes_and_intersections.size());

    for (auto cross_section_i : utility::indices(cs_planes_and_intersections))
    {
        auto& cs = cs_planes_and_intersections[cross_section_i];

        // The (sometimes approximate) intersection of each fiber with the cross-sectional
        // plane has been found.  Now, points need to be placed around the intersections to define
        // a boundary for the cross-section.
        auto intersection_convex_hull_indices =
          simple_math::convex_hull_for_3d_points_on_plane(cs.normal,
                                                          cs.point,
                                                          cs.fiber_intersections);

        std::vector<bool> is_fiber_convex_hull(cs.fiber_intersections.size(), false);
        for (auto& i : intersection_convex_hull_indices)
        {
            is_fiber_convex_hull[i] = true;
        }

        // Transform to 2D
        auto fiber_intersections_2d =
          simple_math::plane_2d_projection(cs.normal, cs.point, cs.fiber_intersections);

        // Fit piecewise spline
        // TODO

        // With the convex hull of each cross-section determined, fit each convex hull with a
        // spline line, which will be used to sample discrete points to define the cross-section
        std::vector<simple_math::Vector3D> convex_hull_points;
        for (auto& i : intersection_convex_hull_indices)
        {
            convex_hull_points.push_back(cs.fiber_intersections[i]);
        }
        convex_hull_points.push_back(cs.fiber_intersections[intersection_convex_hull_indices[0]]);

        // Note: if sparse points result from the convex hull algorithm, then the spline fit can be
        // quite noisy with inflections that are not in the data.  Consequently, interpolated points
        // will be added to fill in large gaps in the points.  We want to interpolate each segment
        // of the convex hull have a roughly equal number of points.
        double min_dist = std::numeric_limits<double>::max();
        for (auto i : utility::range(convex_hull_points.size() - 1))
        {
            auto dist = (convex_hull_points[i] - convex_hull_points[i + 1]).norm();
            if (dist < min_dist)
            {
                min_dist = dist;
            }
        }
        std::vector<simple_math::Vector3D> points_for_spline_fit = { convex_hull_points[0] };
        for (auto i : utility::range(convex_hull_points.size() - 1))
        {
            auto dist = (convex_hull_points[i] - convex_hull_points[i + 1]).norm();
            if (dist > 1.5 * min_dist)
            {
                auto n_segments = int(dist / (1.5 * min_dist));
                auto s          = simple_math::linspace(0.0, 1.0, n_segments);
                for (auto j : utility::range(1, n_segments))
                {
                    auto p = s[j] * (convex_hull_points[i + 1] - convex_hull_points[i]) +
                             convex_hull_points[i];
                    points_for_spline_fit.push_back(p);
                }
            }
            points_for_spline_fit.push_back(convex_hull_points[i + 1]);
        }


        std::vector<simple_math::Point> spline_fit_points;
        auto n_repeated_points = std::min(std::size_t(0), points_for_spline_fit.size() - 1);
        for (auto i : utility::range(points_for_spline_fit.size() - n_repeated_points - 2,
                                     points_for_spline_fit.size() - 1))
        {
            spline_fit_points.push_back(simple_math::to_point(points_for_spline_fit[i]));
        }
        for (auto& p : points_for_spline_fit)
        {
            spline_fit_points.push_back(simple_math::to_point(p));
        }
        for (auto i : utility::range(n_repeated_points))
        {
            spline_fit_points.push_back(simple_math::to_point(points_for_spline_fit[i]));
        }

        // DEBUG TODO delete
        if (cross_section_i == 4176)
        {
            mesh::Mesh m;
            for (auto p : spline_fit_points)
            {
                m.nodes().add(p);
            }
            mesh::write_vtk(m, "spline_points.vtk");
        }

        simple_math::SplineLine3D spline(spline_fit_points,
                                         n_repeated_points,
                                         spline_fit_points.size() - n_repeated_points - 1);

        auto s = simple_math::linspace(0., 1., n_points_around_each_cross_section);

        std::vector<simple_math::Vector3D> discrete_cross_section_points(
          n_points_around_each_cross_section);
        // Note: we leave off s==1 since it should the same point as s==0
        for (auto i : utility::range(s.size() - 1))
        {
            auto sample = spline.at_s(s[i]);
            if (std::isnan(sample.point[0]))
            {
                mesh::Mesh m;
                for (auto p : spline_fit_points)
                {
                    m.nodes().add(p);
                }
                mesh::write_vtk(m, "error_spline_points.vtk");
                _UTIL_ASSERT(!std::isnan(sample.point[0]), "Spline failed");
            }
            discrete_cross_section_points[i] = simple_math::point_as_vector(sample.point);
        }

        discretized_cross_sections[cross_section_i] = { cs.point,
                                                        cs.normal,
                                                        cs.fiber_intersections,
                                                        cs.fiber_intersected_segments,
                                                        cs.fiber_segment_intersection_fractions,
                                                        is_fiber_convex_hull,
                                                        discrete_cross_section_points };
    }

    return discretized_cross_sections;
}

}