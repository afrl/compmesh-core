#pragma once

#include "Common.h"
#include "simple_math/EigenDefines.h"
#include "mesh/Mesh.h"

#include <string>

namespace composite_mesher
{
/// @brief Describes center line of a fiber
struct Fiber
{
    std::vector<simple_math::Point> points;
};

/// @brief Describes center lines of a bundle of fibers
struct FiberBundle
{
    int material_id;
    double diameter;
    std::vector<Fiber> fibers;
};

/// @brief Reads a .fib file that is exported by VTMS
FiberBundle read_fib_file(const std::string& file_path);

/// @brief Convert the fiber bundle to a Mesh of line elements
mesh::Mesh to_mesh(const FiberBundle& bundle);

/// @brief A simple struct to hold the points for a fiber bundle in the form of a matrix which can
/// be convenient for some operations, such as computing winding numbers.
struct PointMatrix
{
    Eigen::MatrixX3d points;
    std::vector<Index> point_to_fiber_id_map;
    std::vector<Index> point_to_index_in_fiber_map;
};

PointMatrix extract_points_matrix(const FiberBundle& bundle);

/// @brief Stores information regarding how far along the fiber each segment is.
struct FiberDistanceInfo
{
    double total_length;
    simple_math::DynamicVector<double> distance_along_line;
    simple_math::DynamicVector<double> fractional_distance_along_line;
};

FiberDistanceInfo get_distance_info(Eigen::Map<const simple_math::DynamicDoubleMatrix> points_ref);

/// @brief Given a segment index and the fraction along that segment, this returns the equivalent
/// distance along the entire fiber.
double compute_distance_along_fiber(const Index& intersected_segment,
                                    const double& segment_intersection_fraction,
                                    const FiberDistanceInfo& distance_info);

/// @brief Given segment indices and fractions along those segments, this returns the equivalent
/// distance along each of the fibers.
std::vector<double> compute_distance_along_fibers(
  const std::vector<Index>& intersected_segments,
  const std::vector<double>& segment_intersection_fractions,
  const std::vector<FiberDistanceInfo>& distance_infos);

/// @brief Given segment indices and fractions along those segments, this returns the equivalent
/// distance divided by the total length of each of the fibers.
std::vector<double> compute_fractions_along_fibers(
  const std::vector<Index>& intersected_segments,
  const std::vector<double>& segment_intersection_fractions,
  const std::vector<FiberDistanceInfo>& distance_infos);

}