#pragma once

#include "Common.h"
#include "FiberBundle.h"
#include "TowSurface.h"
#include "mesh/Mesh.h"
#include "simple_math/EigenDefines.h"

namespace composite_mesher
{
/// @brief Primary function to create a tow surface based on a given bundle of fibers.
///
/// IMPORTANT: the convex hull points need to be enlarged by d/2
///
/// @param bundle - a bundle of fibers that will be used to create the tow surface
/// @param nominal_n_cross_sections - number of cross-sections that would nominally be created if
///        the surface was not adaptively refined in regions of high curvature.
/// @param n_points_around_each_cross_section - number of points to sample around the spline fit of
///        each cross-sectional boundary (determines the resolution of the discretized
///        cross-sections)
/// @param max_angle_allowed_between_fibers_and_plane - maximum angle allowed between any fiber and
///        a vector normal to the plane for a given cross-section (the cross-section will be skipped
///        if the criteria would be violated).
/// @param min_angle_allowed_between_cross_sections - minimum angle allowed between any two
/// sequential cross-sections (more
///        cross-sections will be added if a new cross-section would violate it).
/// @param n_cutbacks_allowed - indicates how many times the increment can be halved when the angle
///        criterion (see parameter above) is violated.
/// @return A triangular mesh of the tow surface (currently without end caps)
VTMSTowSurface generate_surface(const FiberBundle& bundle,
                                const int nominal_n_cross_sections,
                                const int n_points_around_each_cross_section,
                                const double max_angle_allowed_between_fibers_and_plane,
                                const double min_angle_allowed_between_cross_sections,
                                const int n_cutbacks_allowed);

// Below are all of the building blocks for generate_surface, but they are available to use for
// other purposes.

/// @brief See identify_cross_sectional_plane
struct CrossSectionPlane
{
    simple_math::Vector3D point;
    simple_math::Vector3D normal;
    // For each fiber, indicates the index of the segment that is the lowest for which an
    // intersection could be plausible (typically the segment intersected by the last cross-section)
    std::vector<Index> fiber_segments_trust_region_start;
    // For each fiber, indicates the index of the segment that is the highest for which an
    // intersection could be plausible (twice the distance of the extrapolation)
    std::vector<Index> fiber_segments_trust_region_end;
    // For each fiber, the segment index that straddles the length increment for this cross-section.
    // Note: this is also the segment used in the averaging calculation that computes the plane's
    // point and normal.
    std::vector<Index> extrapolated_fiber_intersected_segments;
    std::vector<simple_math::Vector3D> fiber_tangents;
};


/// @brief Identifies the next cross-sectional plane for a given increment the fibers.
///
/// @param bundle - bundle of fibers
/// @param fiber_point_matrices - the points along each fiber, each viewed as an Eigen matrix
/// @param distance_infos - precalculated information about the distance along each fiber for each
///        point in the respective fiber
/// @param delta_s - the relative increment to take along each fiber from the last
///        cross-sectional plane identified
/// @param prev_fiber_intersected_segments - index of the segment along each fiber that was
///        intersected by the last cross-sectional plane identified
/// @param prev_fiber_segment_intersection_fractions - fraction along the segment indicated by
///        prev_fiber_intersected_segments that is the point of intersection with the last
///        cross-sectional plane identified
///
/// @returns A CrossSectionPlane object, which defines the cross-sectional plane with a point and
/// normal and includes the segments in each fiber that are close to the cross-sectional plane
/// based on an extrapolation from the last cross-section identified.
CrossSectionPlane identify_cross_sectional_plane(
  const FiberBundle& bundle,
  const std::vector<Eigen::Map<const simple_math::DynamicDoubleMatrix>>& fiber_point_matrices,
  const std::vector<FiberDistanceInfo>& distance_infos,
  const double& delta_s,
  const std::vector<Index>& prev_fiber_intersected_segments,
  const std::vector<double>& prev_fiber_segment_intersection_fractions);


struct CrossSectionIntersections
{
    std::vector<simple_math::Vector3D> fiber_intersections;
    std::vector<Index> fiber_intersected_segments;
    std::vector<double> fiber_segment_intersection_fractions;
};

/// @brief Finds the intersections between each fiber and the next identified
/// cross-sectional plane.
///
/// @param bundle - bundle of fibers
/// @param fiber_point_matrices - the points along each fiber, each viewed as an Eigen matrix
/// @param plane - point and normal defining next cross-sectional plane and the segments in each
///        fiber that are close to the cross-sectional plane based on an extrapolation from the last
///        cross-section identified, @see identify_cross_sectional_plane.
///
/// @returns An IntersectedCrossSection object, list the coordinate of intersection for each fiber
/// the next cross-sectional plane, the segment along each fiber that is intersected, and the
/// fraction along each of those fibers at which the intersection occurs.
CrossSectionIntersections find_intersections_with_cross_sectional_plane(
  const FiberBundle& bundle,
  const std::vector<Eigen::Map<const simple_math::DynamicDoubleMatrix>>& fiber_point_matrices,
  const CrossSectionPlane& plane);


struct IntersectedCrossSection
{
    simple_math::Vector3D point;
    simple_math::Vector3D normal;
    std::vector<simple_math::Vector3D> fiber_intersections;
    std::vector<Index> fiber_intersected_segments;
    std::vector<double> fiber_segment_intersection_fractions;
    std::vector<double> fiber_intersection_length;
};



/// @brief Iteratively march along the given bundle of fibers and 1) identify plane that
/// will define the next cross-section and 2) find the intersections of that plane with each digital
/// fiber.
///
/// @param bundle - bundle of fibers
/// @param n_cross_sections - number of cross-sections to sample (Notes: #1 must be > 2, #2 can
///        result in a very poor discretization if too few cross-sections are used, and #3 considers
///        extrapolation along each fiber while finding intersections to avoid issues with complex
///        tow paths)
/// @param min_angle_allowed - minimum angle allowed between any two sequential cross-sections (more
///        cross-sections will be added if a new cross-section would violate it).
/// @param n_cutbacks_allowed - indicates how many times the increment can be halved when the angle
///        criterion (see parameter above) is violated.
///
/// @returns A list of cross-sections (defined by a point and normal) and information about the
/// intersection of each fiber with each cross-section.
///
/// Notes for future development:
/// Consider finding bounding box of intersection points, if it changes drastically (% of
/// last one) then back off on increment.
std::vector<IntersectedCrossSection> identify_cross_section_planes_and_intersections(
  const FiberBundle& bundle,
  const int n_cross_sections,
  const double max_angle_deviation_for_cross_section_allowed,
  const double min_angle_between_cross_sections_allowed,
  const int n_cutbacks_allowed);

struct DiscretizedCrossSection
{
    simple_math::Vector3D point;
    simple_math::Vector3D normal;
    std::vector<simple_math::Vector3D> fiber_intersections;
    std::vector<Index> fiber_intersected_segments;
    std::vector<double> fiber_segment_intersection_fractions;
    std::vector<bool> is_fiber_convex_hull;
    std::vector<simple_math::Vector3D> discrete_cross_section_points;
};

/// @brief TODO
///
/// IMPORTANT: the convex hull points need to be enlarged by d/2
///
/// @param cs_planes_and_intersections - a result from calling
///        identify_cross_section_planes_and_intersections
/// @param n_points_around_each_cross_section - number of points to sample around the spline fit of
///        each cross-sectional boundary (determines the resolution of the discretized
///        cross-sections)
///
/// @returns TODO
std::vector<DiscretizedCrossSection> discretize_cross_sections(
  const std::vector<IntersectedCrossSection>& cs_planes_and_intersections,
  const int n_points_around_each_cross_section);


}