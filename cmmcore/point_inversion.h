#ifndef CMMCORE_POINT_INVERSION_H
#define CMMCORE_POINT_INVERSION_H

#include "newthon2.h" // For the Newton's method implementation
#include "nurbs.h" // For NURBS classes
#include <limits>

namespace cmmcore {

/**
 * Helper function to solve a 2x2 linear system
 * 
 * @param matrix 2x2 coefficient matrix
 * @param y right hand side vector
 * @param result solution vector
 * @return true if system was solved successfully
 */
bool solve2x2(const std::array<std::array<double, 2>, 2>& matrix, 
             const std::array<double, 2>& y,
             std::array<double, 2>& result) noexcept {
    double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    if (std::abs(det) < std::numeric_limits<double>::epsilon()) {
        return false;
    }
    // Solve using Cramer's rule
    result[0] = (y[0] * matrix[1][1] - matrix[0][1] * y[1]) / det;
    result[1] = (matrix[0][0] * y[1] - y[0] * matrix[1][0]) / det;
    return true;
}

/**
 * Point inversion for NURBS curves
 * 
 * @param curve NURBS curve to project onto
 * @param P Point to project
 * @param u0 Initial parameter guess
 * @param tol1 Tolerance for Euclidean distance
 * @param tol2 Tolerance for zero cosine measure
 * @param max_iter Maximum number of Newton iterations
 * @return Parameter value that gives closest point on curve
 */
double point_inversion_curve(const NURBSCurve& curve, const vec3& P, 
                           double u0, double tol1, double tol2, 
                           int max_iter = 100) {
    // Temporary storage for evaluated points/derivatives 
    vec3 C_point, C_prime, C_double_prime;
    
    // Get initial point on curve and derivatives
    curve.evaluate(u0, C_point);
    curve.derivative(u0, C_prime);
    curve.second_derivative_fdm(u0, C_double_prime);
    
    double ui = u0;
    
    for (int i = 0; i < max_iter; ++i) {
        // Evaluate f(u) = C'(u) · (C(u) - P)
        double fi = C_prime.dot(C_point - P);
        
        // Evaluate f'(u) = C''(u) · (C(u) - P) + C'(u) · C'(u)
        double fpi = C_double_prime.dot(C_point - P) + C_prime.dot(C_prime);
        
        // Check if derivative is too small
        if (std::abs(fpi) < std::numeric_limits<double>::epsilon()) {
            break;
        }
        
        // Newton step
        double ui1 = ui - fi / fpi;
        
        // Convergence check 1: Parameter change small enough
        if (std::abs(ui1 - ui) * C_prime.length() <= tol1) {
            break;
        }
        
        // Update point and derivatives at new parameter
        curve.evaluate(ui1, C_point);
        
        // Convergence check 2: Distance to point small enough 
        if ((C_point - P).length() <= tol1) {
            break;
        }
        
        curve.derivative(ui1, C_prime);
        
        // Convergence check 3: Angle between difference and tangent small enough
        double cos_theta = std::abs(C_prime.dot(C_point - P)) / 
                          (C_prime.length() * (C_point - P).length());
        if (cos_theta <= tol2) {
            break;
        }
        
        curve.second_derivative_fdm(ui1, C_double_prime);
        ui = ui1;
    }
    
    return ui;
}

/**
 * Point inversion for NURBS surfaces
 * 
 * @param surface NURBS surface to project onto
 * @param P Point to project
 * @param u0 Initial u parameter guess
 * @param v0 Initial v parameter guess 
 * @param tol1 Tolerance for Euclidean distance
 * @param tol2 Tolerance for zero cosine measure
 * @param max_iter Maximum number of Newton iterations
 * @return Pair of (u,v) parameters that give closest point on surface
 */
std::array<double,2> point_inversion_surface(const NURBSSurface& surface, const vec3& P,
                                           double u0, double v0, 
                                           double tol1, double tol2,
                                           int max_iter = 100) {
    // Initialize current parameter point
    std::array<double,2> uivi = {u0, v0};
    std::array<double,2> delta = {0.0, 0.0};
    
    // Temporary storage for evaluated points/derivatives
    vec3 S_point, du, dv;
    
    for (int i = 0; i < max_iter; ++i) {
        // Get point and derivatives at current parameters 
        surface.evaluate(uivi[0], uivi[1], S_point);
        surface.derivative_u(uivi[0], uivi[1], du);
        surface.derivative_v(uivi[0], uivi[1], dv);
        
        // Computing f(u,v) = du · (S(u,v) - P)
        double f_val = du.dot(S_point - P);
        // Computing g(u,v) = dv · (S(u,v) - P) 
        double g_val = dv.dot(S_point - P);
        
        // Set up Jacobian matrix
        std::array<std::array<double,2>,2> J = {{
            {du.dot(du), du.dot(dv)},
            {dv.dot(du), dv.dot(dv)}
        }};
        
        // Right hand side vector
        std::array<double,2> k = {f_val, g_val};
        
        // Solve linear system J * delta = -k
        bool success = solve2x2(J, k, delta);
        if (!success) {
            throw std::runtime_error("Failed to solve 2x2 system in surface inversion");
        }
        
        // Update parameters
        std::array<double,2> ui1vi1 = {
            uivi[0] - delta[0],
            uivi[1] - delta[1]
        };
        
        // Check convergence criteria
        
        // 1. Parameter update small enough?
        if ((ui1vi1[0] - uivi[0]) * du.length() + 
            (ui1vi1[1] - uivi[1]) * dv.length() <= tol1) {
            break;
        }
        
        // 2. Distance to point small enough?
        surface.evaluate(ui1vi1[0], ui1vi1[1], S_point);
        if ((S_point - P).length() <= tol1) {
            break;
        }
        
        // 3. Perpendicularity condition met?
        surface.derivative_u(ui1vi1[0], ui1vi1[1], du);
        vec3 pt_diff = S_point - P;
        if (std::abs(du.dot(pt_diff)) / 
            (du.length() * pt_diff.length()) <= tol2) {
            break;
        }
        
        uivi = ui1vi1;
    }
    
    return uivi;
}

}  // namespace cmmcore

#endif // CMMCORE_POINT_INVERSION_H