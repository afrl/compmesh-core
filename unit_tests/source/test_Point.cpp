#include "simple_math/Utility.h"

#include <iostream>

using namespace simple_math::point_operators;

int main()
{
    simple_math::Point a = { 0, 1, 2 };
    simple_math::Point b = { 2, 3, 4 };

    auto c = a + b;
    assert(c == simple_math::Point({ 2, 4, 6 }));

    auto d = 2.0 * a;
    assert(d == simple_math::Point({ 0, 2, 4 }));

    return 0;
}