#include "simple_math/ConvexHull.h"

#include <iostream>

using namespace simple_math;

int main()
{
    // ^ y
    // |
    // p3        p1
    // |
    // |    p4
    // |  X <- plane point
    // p0 ------ p2  -- x
    std::vector<Vector3D> points = { { 0, 0, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0.5, 0.5, 0 } };

    auto convex_hull_point_indices =
      convex_hull_for_3d_points_on_plane({ 0, 0, 1 }, { 0.25, 0.25, 0. }, points);
    assert(convex_hull_point_indices == std::vector<std::size_t>({ 0, 2, 1, 3 }));

    convex_hull_point_indices =
      convex_hull_for_3d_points_on_plane({ 0, 0, -1 }, { 0.25, 0.25, 0. }, points);
    assert(convex_hull_point_indices == std::vector<std::size_t>({ 0, 3, 1, 2 }));

    return 0;
}