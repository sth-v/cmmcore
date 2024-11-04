#ifndef CMMCORE_BOUNDARY_INTERSECTION_H
#define CMMCORE_BOUNDARY_INTERSECTION_H

#include "nurbs.h" // For NURBS classes
#include <vector>
#include <array>
#include <optional>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <limits>
#include <string>
#include <set>
#include "csx.h" // For nurbs_csx

namespace cmmcore {

// Forward declarations
class NURBSCurve;
class NURBSSurface;

/**
 * Extract an isocurve from a NURBS surface at a given parameter in the u or v direction
 */
NURBSCurve extract_isocurve(const NURBSSurface& surface, double param, const std::string& direction = "u") {

    if (direction != "u" && direction != "v") {
        throw std::invalid_argument("Direction must be either 'u' or 'v'.");
    }

    auto interval = surface.interval();
    
    if (direction == "u") {
        // For u-direction: we fix u and vary v
        auto param_range = interval[0];  // u range
        if (param < param_range[0] || param > param_range[1]) {
            throw std::invalid_argument("Parameter " + std::to_string(param) + 
                                      " is out of range [" + std::to_string(param_range[0]) + 
                                      "," + std::to_string(param_range[1]) + "]");
        }

        // Find span and basis functions in u direction
        int n_u = surface._size[0] - 1;
        int degree_u = surface._degree[0];
        int span = find_span(n_u, degree_u, param,  surface._knots_u, false);
        std::array<double,CMMCORE_DEG_STACK_SIZE> basis{};
        basis_funs(span, param, degree_u, surface._knots_u, basis);


        // The curve will have as many control points as the surface has in v direction
        size_t m = surface._size[1];
        NURBSCurve crv(m, surface._degree[1]);
        for (size_t i=0;i<surface._knots_v.size();i++) {
            crv.knots[i]=surface._knots_v[i];
        }



        // Compute control points for the extracted curve
        for (int i = 0; i < m; ++i) {
            vec4 point(0.0);
            for (int j = 0; j <= degree_u; ++j) {
                point +=  surface._control_points[span - degree_u + j][i] * basis[j];
            }
            crv.control_points[i].set(point);
        }

        return crv;
    } 
    else {  // direction == "v"
        // For v-direction: we fix v and vary u
        auto param_range = interval[1];  // v range
        if (param < param_range[0] || param > param_range[1]) {
            throw std::invalid_argument("Parameter " + std::to_string(param) + 
                                      " is out of range [" + std::to_string(param_range[0]) + 
                                      "," + std::to_string(param_range[1]) + "]");
        }

        // Find span and basis functions in v direction
        int n_v = surface._size[1] - 1;
        int degree_v = surface._degree[1];
        int span = find_span(n_v, degree_v, param, surface._knots_v, false);
        std::array<double,CMMCORE_DEG_STACK_SIZE> basis{};
        basis_funs(span, param, degree_v, surface._knots_v,basis);


        // The curve will have as many control points as the surface has in u direction
        int m = surface._size[0];
        NURBSCurve crv(m, surface._degree[0]);

        // Compute control points for the extracted curve
        for (int i = 0; i < m; ++i) {
            vec4 point(0.0);
            for (int j = 0; j <= degree_v; ++j) {
                point += surface._control_points[i][span - degree_v + j] * basis[j];
            }
            crv.control_points[i] = point;
        }

        return crv;
    }
}

/**
 * Extract the four boundary curves of a NURBS surface
 */
inline void extract_surface_boundaries(const NURBSSurface& surface, std::array<NURBSCurve,4> &result) {
    auto& interval = surface._interval;
    double u_min = interval[0][0], u_max = interval[0][1];
    double v_min = interval[1][0], v_max = interval[1][1];
    
    // Extract iso-curves at the boundaries
    result[0] = std::move(extract_isocurve(surface, u_min, "u"));  // v-direction curve at u=0
      result[1] = std::move(extract_isocurve(surface, u_max, "u"));  // v-direction curve at u=1
     result[2]= std::move(extract_isocurve(surface, v_min, "v"));  // u-direction curve at v=0
     result[3] = std::move(extract_isocurve(surface, v_max, "v"));  // u-direction curve at v=1
    

}

/**
 * Class representing an intersection point between a boundary curve and a surface
 */
struct IntersectionPoint {

    vec3 point;
    double curve_param;
    std::array<double, 2> surface_params;
    std::array<double, 2> surface1_params;
    std::array<double, 2> surface2_params;
    int boundary_index;
    bool is_from_first_surface;
    std::array<std::array<double, 2>, 2> interval;
    IntersectionPoint(const vec3& point, 
                     double curve_param,
                     const std::array<double, 2>& surface_params,
                     int boundary_index,
                     bool is_from_first_surface,
                     const std::array<std::array<double, 2>, 2>& interval)
        : point(point)
        , curve_param(curve_param)
        , surface_params(surface_params)
        , boundary_index(boundary_index)
        , is_from_first_surface(is_from_first_surface)
        , interval(interval) {
        
        // Store parameters for both surfaces
        if (is_from_first_surface) {
            surface1_params = boundary_index_to_params(boundary_index, curve_param);
            surface2_params = surface_params;
        } else {
            surface1_params = surface_params;
            surface2_params = boundary_index_to_params(boundary_index, curve_param);
        }
    }

    std::pair<std::array<double, 2>, std::array<double, 2>> get_start_params() const {
        if (is_from_first_surface) {
            // For first surface, convert boundary index to fixed parameter
            auto [u1, v1] = boundary_index_to_params(boundary_index, curve_param);
            // For second surface, use the found parameters
            double u2 = surface_params[0], v2 = surface_params[1];
            return {{u1, v1}, {u2, v2}};
        } else {
            // For second surface, convert boundary index to fixed parameter
            auto [u2, v2] = boundary_index_to_params(boundary_index, curve_param);
            // For first surface, use the found parameters
            double u1 = surface_params[0], v1 = surface_params[1];
            return {{u1, v1}, {u2, v2}};
        }
    }

    //const vec3& point() const { return point_; }
    //double curve_param() const { return curve_param_; }
    //const std::array<double, 2>& surface_params() const { return surface_params_; }
    //int boundary_index() const { return boundary_index_; }
    //bool is_from_first_surface() const { return is_from_first_surface_; }

private:
     std::array<double, 2> boundary_index_to_params(int boundary_index, double param) const {
        if (boundary_index == 0)      return {interval[0][0], param}; // u=0 curve
        else if (boundary_index == 1) return {interval[0][1], param}; // u=1 curve
        else if (boundary_index == 2) return {param, interval[1][0]}; // v=0 curve
        else                         return {param, interval[1][1]}; // v=1 curve
    }


};


/**
 * Find all intersection points between the boundaries of two NURBS surfaces
 */
inline void find_boundary_intersections(
    const NURBSSurface& surf1, 
    const NURBSSurface& surf2,
    std::vector<IntersectionPoint>& intersection_points,
    double tol = 1e-6) {
    

    
    // Get boundaries of both surfaces
    std::array<NURBSCurve,4> boundaries1;
    extract_surface_boundaries(surf1,boundaries1);
    std::array<NURBSCurve,4> boundaries2;
    extract_surface_boundaries(surf2,boundaries2);
    CSXIntersections intersections{};
    // Find intersections of surf1's boundaries with surf2
    for (size_t i = 0; i < boundaries1.size(); ++i) {
        intersections.clear();
        csx(boundaries1[i], surf2, intersections,tol);

        for (const auto& [intersection_type, point, params] : intersections) {
            intersection_points.emplace_back(
                point,
                params[0],  // curve parameter
                std::array<double,2>{params[1], params[2]}, // surface parameters
                i,
                true,
                surf1._interval
            );
        }
    }
    
    // Find intersections of surf2's boundaries with surf1
    for (size_t i = 0; i < boundaries2.size(); ++i) {
        intersections.clear();
        csx(boundaries2[i], surf1, intersections,tol);
        for (const auto& [intersection_type, point, params] : intersections) {
            intersection_points.emplace_back(
                point,
                params[0],  // curve parameter
                std::array<double,2>{params[1], params[2]}, // surface parameters
                i,
                false,
                surf2._interval
            );
        }
    }
    
    // Remove duplicate points (within tolerance)

    auto unique_points=Vec3Set();
    const size_t sx=intersection_points.size();
    for (int i = sx- 1; i >= 0; --i) {


        auto& point=intersection_points[i];

        if (!contains(unique_points, point.point))
        {
            unique_points.emplace(point.point);
            continue;

        }

        intersection_points.erase(intersection_points.begin()+i);


    }
    
    //return unique_points;
}

/**
 * Sort boundary intersection points into connected sequences that form intersection curves
 */
inline void sort_boundary_intersections(
    const std::vector<IntersectionPoint>& points,std::vector<std::vector<IntersectionPoint>>& result) {
    
    if (points.empty()) {
        return ;
    }
    
    // Start with all points in unassigned set
    std::set<size_t> unassigned;
    for (size_t i = 0; i < points.size(); ++i) {
        unassigned.insert(i);
    }
    
    //std::vector<std::vector<IntersectionPoint>> sequences;
    
    while (!unassigned.empty()) {
        // Start a new sequence with the first unassigned point
        std::vector<IntersectionPoint> current_sequence;
        size_t start_idx = *unassigned.begin();
        unassigned.erase(start_idx);
        current_sequence.push_back(points[start_idx]);
        
        // Try to find the next closest point
        while (true) {
            const auto& current_point = current_sequence.back();
            std::optional<size_t> nearest_idx;
            double min_distance = std::numeric_limits<double>::infinity();
            
            // Find the closest unassigned point
            for (size_t idx : unassigned) {
                const auto& candidate = points[idx];
                double distance = (current_point.point - candidate.point).length();
                
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_idx = idx;
                }
            }
            
            // If no close point found or sequence has 2 points, end sequence
            if (!nearest_idx || current_sequence.size() == 2) {
                break;
            }
            
            // Add the nearest point to sequence and remove from unassigned
            current_sequence.push_back(points[*nearest_idx]);
            unassigned.erase(*nearest_idx);
        }
        
        result.push_back(std::move(current_sequence));
    }
    

}

} // namespace cmmcore

#endif // CMMCORE_BOUNDARY_INTERSECTION_H