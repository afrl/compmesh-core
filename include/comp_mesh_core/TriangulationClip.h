#pragma once

#include "Triangulation.h"

namespace composite_mesher
{

/// @brief Clips a triangulation using a plane (defined by a point and normal vector), keeping the
///        surface that lies on the positive side of the plane.
void clip_by_plane(Triangulation& t,
                   const simple_math::Point& plane_point,
                   const simple_math::Point& plane_normal);

}