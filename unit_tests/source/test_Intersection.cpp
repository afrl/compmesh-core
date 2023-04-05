#include "Eigen/src/Core/util/DisableStupidWarnings.h"
#include "simple_math/Intersection.h"
#include "simple_math/Utility.h"
#include "unsupported/Eigen/NonLinearOptimization"

#include <iomanip>
#include <iostream>

static constexpr double tol        = 1e-6;
static constexpr bool print_status = false;

/// @brief A functor for a paraboloid to use with non-linear solvers to find the signed distance and
/// nearest point to a given point.
///
/// Defines a paraboloid of the form z(x, y) = x^2 + y^2.  Parametric coordinates are used (phi,
/// psi) transforming the form to: x(phi, psi) = phi, y(phi, psi) = psi, and z(phi, psi) = phi^2 +
/// psi^2.
///
/// Adopted from the example in eigen/unsupported/test/NonLinearOptimization.cpp
struct ParaboloidSignedDistanceFunctor
{
    using Scalar       = double;
    using InputType    = Eigen::VectorXd;
    using ValueType    = Eigen::VectorXd;
    using JacobianType = Eigen::MatrixXd;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    static const int n_inputs_ = 3;
    static const int n_values_ = 3;
    const Eigen::VectorXd point_q_;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Methods

    /// @param q - the point at which to find the closest point on the surface.
    ParaboloidSignedDistanceFunctor(const Eigen::VectorXd& q)
      : point_q_(q)
    {
    }

    /// @param inputs:
    /// { psi, phi, s } where psi and phi are parametric coordinates for the paraboloid and s if the
    /// signed distance between the given point, q, and the surface.
    ///
    /// @param values:
    /// equivalent to the residual vector, R, for the system of equations:
    /// R = p + s*n - q,
    /// where q is the given point for which to find the closest point on the surface, p, s is the
    /// signed distnace between p and q, and n is the normal to the surface at p.
    ///
    /// R should be { 0, 0, 0 } at any point on the surface for which a line passes through
    /// p and q parallel to the normal.  For a paraboloid, this should be unique due to the
    /// constant sign of the curvature (remains concave for all values of psi and phi).
    int operator()(const Eigen::VectorXd& inputs, Eigen::VectorXd& values) const
    {
        assert(inputs.size() == 3);
        assert(values.size() == 3);

        const double& phi = inputs[0];
        const double& psi = inputs[1];
        const double& s   = inputs[2]; // signed distance

        auto p = eval_p(phi, psi);
        auto n = eval_n(phi, psi);
        values = p + (s * n) - point_q_;

        if (print_status)
        {
            const auto w = 12;
            std::cout << "   phi = " << std::setw(w) << phi << "   psi = " << std::setw(w) << psi
                      << "   s = " << std::setw(w) << s << "   |R| = " << std::setw(w)
                      << values.norm() << std::endl;
        }

        return 0;
    }

    int inputs() const
    {
        return n_inputs_;
    }
    int values() const
    {
        return n_values_;
    }

    /// @brief Evaluates the global coordinate of a point on the surface given the parametric point
    /// { phi, psi }.
    static Eigen::Vector3d eval_p(const double& phi, const double& psi)
    {
        return { phi, psi, phi * phi + psi * psi };
    }

    /// @brief Evaluates the unit vector that is normal to the paraboloid at the parametric point
    /// { phi, psi }.
    ///
    /// Let t_phi be the vector tangent to the surface along the phi axis, which is defined by
    /// \nabla p(phi, psi), and let t_psi be the vector tangent to the surface along the psi axis,
    /// which is defined by \nabla p(phi, psi).  The unit normal n would then be given by:
    /// n = ( t_phi x t_psi ) / sqrt( |t_phi|^2 * |t_psi|^2 + |t_phi \dot t_psi|^2 )
    /// However, this is rearranged a bit to make the evaluation a bit faster in the implementation.
    static Eigen::Vector3d eval_n(const double& phi, const double& psi)
    {
        const auto n_den = sqrt(1 + 4 * phi * phi + 4 * psi * psi);
        return { -2 * phi / n_den, -2 * psi / n_den, 1 / n_den };
    }

    /// @brief A helper for testing if a point is actually on the surface.
    static bool is_point_on_surface(simple_math::Point p)
    {
        // Fortunately, there is a trivial inverse map for (x and phi) and (y and psi)
        const double phi = p[0];
        const double psi = p[1];
        const double z   = phi * phi + psi * psi;
        return abs(z - p[2]) < tol;
    }
};

struct ParaboloidNearestPoint
{
    /// @brief Nearest point on the surface in the global coordinate system.
    Eigen::Vector3d p;
    /// @brief Unit normal vector at nearest point on surface, p.
    Eigen::Vector3d n;
    /// @brief Signed distance (in global coordinate system).
    double s;
};

ParaboloidNearestPoint find_nearest_point(Eigen::Vector3d q)
{
    ParaboloidSignedDistanceFunctor functor(q);
    Eigen::HybridNonLinearSolver<ParaboloidSignedDistanceFunctor> solver(functor);

    Eigen::VectorXd vars(3);
    vars.setConstant(0.0); // initial guess is { 0, 0, 0 }

    if (print_status)
    {
        std::cout << "Solving for point nearest: { " << q[0] << ", " << q[1] << ", " << q[2] << " }"
                  << std::endl;
    }

    auto return_code = solver.hybrd1(vars);
    assert(return_code == 1);
    if (print_status)
    {
        std::cout << std::endl;
    }

    const auto& phi = vars[0];
    const auto& psi = vars[1];
    const auto& s   = vars[2];

    return { ParaboloidSignedDistanceFunctor::eval_p(phi, psi),
             ParaboloidSignedDistanceFunctor::eval_n(phi, psi),
             s };
}

void test_paraboloid_solver(Eigen::Vector3d q, ParaboloidNearestPoint reference_sol)
{
    auto nearest_point = find_nearest_point(q);

    assert((nearest_point.p - reference_sol.p).norm() < tol);
    assert((nearest_point.n - reference_sol.n).norm() < tol);
    assert(abs(nearest_point.s - reference_sol.s) < tol);
}

int main()
{
    // Checks to make sure the paraboloid nearest-point solver is implemented properly

    test_paraboloid_solver({ 0, 0, 0 }, { { 0, 0, 0 }, { 0, 0, 1 }, 0 });
    test_paraboloid_solver({ 0, 0, 2 }, { { 0, 0, 0 }, { 0, 0, 1 }, 2 });
    test_paraboloid_solver({ 1, 1, 2 }, { { 1, 1, 2 }, { -2. / 3., -2. / 3., 1. / 3. }, 0 });
    // Based on answer from a Python script
    test_paraboloid_solver({ 2, 2, 2 },
                           { { 1.09791167, 1.09791167, 2.41082008 },
                             { -0.67306909, -0.67306909, 0.30652242 },
                             -1.3402611090231529 });


    // Test the intersection algorithm

    simple_math::SignedDistanceFunction paraboloid_level_set = {
        [](const simple_math::Point& p) -> double {
            auto nearest_point = find_nearest_point(simple_math::point_as_vector(p));
            return nearest_point.s;
        },
        2
    };

    // Case: single intersection
    {

        simple_math::Point segment_p1 = { 0, 0, 1 };
        simple_math::Point segment_p2 = { 2, 2, 1 };
        auto intersections = simple_math::find_segment_level_set_intersection(paraboloid_level_set,
                                                                              segment_p1,
                                                                              segment_p2);
        assert(intersections.size() == 1);
        assert(intersections[0].type == simple_math::Intersection::Type::INTERIOR);
        assert(ParaboloidSignedDistanceFunctor::is_point_on_surface(intersections[0].point));
    }

    // Case: multiple intersections
    {

        simple_math::Point segment_p1 = { -2, -2, 1 };
        simple_math::Point segment_p2 = { 2, 2, 1 };
        auto intersections = simple_math::find_segment_level_set_intersection(paraboloid_level_set,
                                                                              segment_p1,
                                                                              segment_p2);
        assert(intersections.size() == 2);
        assert(intersections[0].type == simple_math::Intersection::Type::INTERIOR);
        assert(ParaboloidSignedDistanceFunctor::is_point_on_surface(intersections[0].point));
        assert(intersections[1].type == simple_math::Intersection::Type::INTERIOR);
        assert(ParaboloidSignedDistanceFunctor::is_point_on_surface(intersections[1].point));
    }

    // Case: no intersections
    {

        simple_math::Point segment_p1 = { -1, -1, 3 };
        simple_math::Point segment_p2 = { 1, 1, 3 };
        auto intersections = simple_math::find_segment_level_set_intersection(paraboloid_level_set,
                                                                              segment_p1,
                                                                              segment_p2);
        assert(intersections.size() == 0);
    }

    // Case: intersection at tangent point
    {

        simple_math::Point segment_p1 = { 0.5, 0.5, 1 };
        simple_math::Point segment_p2 = { 2, 2, 4 };
        auto intersections = simple_math::find_segment_level_set_intersection(paraboloid_level_set,
                                                                              segment_p1,
                                                                              segment_p2);
        assert(intersections.size() == 1);
        assert(intersections[0].type == simple_math::Intersection::Type::INTERIOR);
        assert(ParaboloidSignedDistanceFunctor::is_point_on_surface(intersections[0].point));
    }

    // Case: intersection at end point
    {

        simple_math::Point segment_p1 = { 1, 1, 2 };
        simple_math::Point segment_p2 = { 2, 2, 4 };
        auto intersections = simple_math::find_segment_level_set_intersection(paraboloid_level_set,
                                                                              segment_p1,
                                                                              segment_p2);
        assert(intersections.size() == 1);
        assert(intersections[0].type == simple_math::Intersection::Type::ENDPOINT);
        assert(ParaboloidSignedDistanceFunctor::is_point_on_surface(intersections[0].point));
    }

    return 0;
}