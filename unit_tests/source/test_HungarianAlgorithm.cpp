#include "simple_math/HungarianAlgorithm.h"

#include <assert.h>

using namespace simple_math;

int main()
{
    // ^ y
    // |
    // p3        p1
    // |
    // |    p4
    // |
    // p0 ------ p2  -- x
    std::vector<Point> point_set_1 = { { 0, 0, 0 },
                                       { 1, 1, 0 },
                                       { 1, 0, 0 },
                                       { 0, 1, 0 },
                                       { 0.5, 0.5, 0 } };

    // ^ y
    // |
    // p4        p2
    // |
    // |    p0
    // |
    // p1 ------ p3  -- x
    std::vector<Point> point_set_2 = { { 0.6, 0.55, 1.1 },
                                       { 0.2, 0, 1.05 },
                                       { 1.1, 0.96, 0.9 },
                                       { 1.05, 0.05, 1.0 },
                                       { 0.1, 1.1, 0.95 } };

    auto assignments = solve_linear_sum_assignment(point_set_1, point_set_2);
    assert(assignments[0] == 1);
    assert(assignments[1] == 2);
    assert(assignments[2] == 3);
    assert(assignments[3] == 4);
    assert(assignments[4] == 0);

    return 0;
}