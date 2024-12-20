//
// Created by Andrew Astakhov on 28.09.24.
//

#ifndef NURBS_H
#define NURBS_H
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>

#include <unordered_set>
#include <fstream>
#include <sstream>

#include "nurbs_utils.h"
#include "utils.h"
#include "aabb.h"

#include "cmmcore.h"
#include "integrate.h"

namespace cmmcore {
#define CMMCORE_NURBS_INTEGRATION_TOL 1e-3

    namespace nurbs_helpers {
        /**
         * Flips the control points grid from U-major to V-major.
         *
         * @param cpts Input control points grid (U-major).
         * @param result Output control points grid (V-major).
         */
        inline void flipControlPointsU(const std::vector<std::vector<vec4> > &cpts,
                                       std::vector<std::vector<vec4> > &result) {
            size_t size_u = cpts.size();
            size_t size_v = cpts[0].size();
            result.clear();


            for (size_t i = 0; i < size_v; i++) {
                std::vector<vec4> new_row{};
                for (size_t j = 0; j < size_u; j++) {
                    new_row.push_back(cpts[j][i]);
                }
                result.push_back(new_row);
            }
        }
    }

    class NURBSCurve {
    public:
        NURBSCurve() = default;

        NURBSCurve(const std::vector<vec4> &control_points): NURBSCurve(
            control_points, control_points.size() >= 4 ? 3 : 1, false) {
        };

        NURBSCurve(const size_t cpt_size, const int degree, const bool periodic = false): _degree(degree),
            _periodic(periodic) {
            control_points.resize(cpt_size);
            if (_periodic) {
                generate_knots_periodic();
            } else {
                generate_knots();
            }

            update_interval();
        };

        NURBSCurve(const std::vector<vec4> &control_points, const bool periodic): NURBSCurve(
            control_points, control_points.size() >= 4 ? 3 : 1, periodic) {
        };

        NURBSCurve(const std::vector<vec3> &control_points, const int degree): NURBSCurve(
            control_points, degree, false) {
        };

        NURBSCurve(const std::vector<vec4> &control_points, const int degree): NURBSCurve(
            control_points, degree, false) {
        };

        // Constructors
        NURBSCurve(const std::vector<vec3> &_control_points, const int degree,
                   const bool periodic): _degree(degree), _periodic(periodic) {
            // Initialize control points
            const size_t nnn = _control_points.size();
            control_points.resize(nnn);
            for (size_t i = 0; i < nnn; ++i) {
                control_points[i].set(_control_points[i]);
            }

            // Initialize knots
            if (_periodic) {
                generate_knots_periodic();
            } else {
                generate_knots();
            }
            if (_periodic) {
                make_periodic();
            }
            rebuildAABB();
        }

        NURBSCurve(const std::vector<vec4> &control_points, const int degree,
                   bool periodic) : control_points(control_points), _degree(degree), _periodic(periodic) {
            // Initialize knots
            if (_periodic) {
                generate_knots_periodic();
            } else {
                generate_knots();
            }
            if (_periodic) {
                make_periodic();
            }
            rebuildAABB();
        }

        NURBSCurve(const std::vector<vec4> &control_points, const int degree,
                   const std::vector<double> &knots, bool periodic)
            : control_points(control_points), _degree(degree), knots(knots), _periodic(periodic) {
            update_interval();
            rebuildAABB();
        }

        // Copy constructor and assignment operator
        NURBSCurve(const NURBSCurve &other) {
            control_points = other.control_points;
            _degree = other._degree;
            knots = other.knots;
            _aabb = other._aabb;
            _periodic = other._periodic;
            _interval = other._interval;
            _hull = other._hull;
        };

        NURBSCurve(NURBSCurve &&other) {
            if (this != &other) {
                control_points = other.control_points;
                _degree = other._degree;
                knots = other.knots;
                _aabb = other._aabb;
                _interval = other._interval;
                _hull = other._hull;
            }
        }

        NURBSCurve &operator=(const NURBSCurve &other) {
            control_points = other.control_points;
            _degree = other._degree;
            knots = other.knots;
            _aabb = other._aabb;
            _periodic = other._periodic;
            _interval = other._interval;
            _hull = other._hull;
            return *this;
        };

        // Methods

        // Check if the NURBS curve is periodic
        constexpr bool is_periodic() const {
            for (int j = 0; j < _degree; ++j) {
                for (int i = 0; i < 4; ++i) {
                    if (control_points[j][i] != control_points[control_points.size() - _degree + j][i]) {
                        return false;
                    }
                }
            }
            return true;
        }

        // Update interval based on knots
        void update_interval() {
            _interval[0] = std::min_element(knots.begin(), knots.end())[0];
            _interval[1] = std::max_element(knots.begin(), knots.end())[0];
        }

        // Hook function called after knots are updated
        void knots_update_hook() {
            update_interval();
        }

        // Normalize the knot vector
        void normalize_knots() {
            double start = std::min_element(knots.begin(), knots.end())[0];
            double end = std::max_element(knots.begin(), knots.end())[0];
            double d = 1.0 / (end - start);

            for (auto &knot: knots) {
                knot = (knot - start) * d;
            }
            knots_update_hook();
        }

        // Generate default open knot vector
        void generate_knots() {
            int n = static_cast<int>(control_points.size());
            knots.clear();
            knots.resize(n + _degree + 1);

            std::fill(knots.begin(), knots.begin() + _degree + 1, 0.0);
            std::iota(knots.begin() + _degree + 1, knots.end() - _degree, 1.0);
            std::fill(knots.end() - _degree, knots.end(), static_cast<double>(n - _degree));

            knots_update_hook();
        }

        // Make the NURBS curve periodic
        void make_periodic() {
            if (is_periodic()) {
                return;
            }
            int n = static_cast<int>(control_points.size());
            int new_n = n + _degree;
            std::vector<vec4> new_control_points(new_n);

            // Copy original control points
            std::copy(control_points.begin(), control_points.end(), new_control_points.begin());

            // Add the first degree control points to the end
            for (int i = 0; i < _degree; ++i) {
                new_control_points[n + i] = control_points[i];
            }
            control_points = std::move(new_control_points);
            generate_knots_periodic();
            _periodic = true;
        }

        // Generate knot vector for a periodic NURBS curve
        void generate_knots_periodic() {
            int n = static_cast<int>(control_points.size());
            int m = n + _degree + 1;
            knots.resize(m);
            for (int i = 0; i < m; ++i) {
                knots[i] = static_cast<double>(i - _degree);
            }
            knots_update_hook();
        }

        // Insert knots into the curve
        void insert_knot(double t, int count) {
            int n = static_cast<int>(control_points.size());
            int new_count = n + count;
            const std::vector<vec4> cpts = control_points;

            // Find knot span
            int span = find_span(n - 1, _degree, t, knots, false);


            // Compute new knot vector
            std::vector<double> k_v = knot_insertion_kv(knots, t, span, count);
            int s_u = find_multiplicity(t, knots);
            // Compute new control points
            control_points.resize(new_count);

            for (int i = 0; i <= count; ++i) {
                control_points[i] = cpts.back() - i;
            }
            knot_insertion(_degree, knots, cpts, t, count, s_u, span, control_points);
            // Update curve
            knots = std::move(k_v);
            update_interval();
        }

        // Split the curve at a given parameter value
        std::pair<NURBSCurve, NURBSCurve> split(const double param, const bool normalize_knots = false) const {
            if (param <= _interval[0] || param >= _interval[1] ||
                std::fabs(param - _interval[0]) <= 1e-12 || std::fabs(param - _interval[1]) <= 1e-12) {
                throw std::invalid_argument("Cannot split from the domain edge.");
            }
            int n_ctrlpts = static_cast<int>(control_points.size());
            auto ks = find_span(n_ctrlpts, _degree, param, knots, false) - _degree + 1;
            int s = find_multiplicity(param, knots);
            int r = _degree - s;
            // Insert knot
            NURBSCurve temp_obj(this->control_points, this->_degree, this->knots, this->_periodic);
            temp_obj.update_interval();
            temp_obj.insert_knot(param, r);
            // Knot span index
            int knot_span = find_span(temp_obj.control_points.size(), _degree, param, temp_obj.knots, false) + 1;

            // Create knot vectors for the two new curves
            std::vector<double> surf1_kv(knot_span);
            // printf("%d\n",knot_span);
            std::vector<double> surf2_kv(temp_obj.knots.size() - knot_span);
            //surf1_kv.assign(temp_knots.begin(), temp_knots.begin()+ knot_span);
            //surf2_kv.assign(temp_knots.end()-knot_span, temp_knots.end());
            for (int i = 0; i < knot_span; ++i) {
                surf1_kv[i] = temp_obj.knots[i];
            }
            for (int i = 0; i < temp_obj.knots.size() - knot_span; ++i) {
                surf2_kv[i] = temp_obj.knots[i + knot_span];
            }
            // Add param to the end of surf1_kv and beginning of surf2_kv
            surf1_kv.push_back(param);
            //surf2_kv.insert(surf2_kv.begin(),param);
            for (int j = 0; j <= _degree; ++j) {
                surf2_kv.insert(surf2_kv.begin(), param);
            }
            // Control points for the new curves
            std::vector<vec4> curve1_ctrlpts(temp_obj.control_points.begin(),
                                             temp_obj.control_points.begin() + ks + r);
            std::vector<vec4> curve2_ctrlpts(temp_obj.control_points.begin() + ks + r - 1,
                                             temp_obj.control_points.end());
            // Create new NURBS curves
            NURBSCurve curve1(curve1_ctrlpts, _degree, surf1_kv, false);
            NURBSCurve curve2(curve2_ctrlpts, _degree, surf2_kv, false);
            if (normalize_knots) {
                curve1.normalize_knots();
                curve2.normalize_knots();
            }

            return {curve1, curve2};
        }

        // Evaluate a point on the NURBS curve at parameter t
        void evaluate(const double t, vec3 &result) const {
            int n = static_cast<int>(control_points.size()) - 1;
            result.set(0, 0, 0);

            // Assume curve_point is already implemented
            curve_point(n, _degree, knots, control_points, t, result, _periodic);
        }

        void derivative(const double t, vec3 &result) const {
            int n = static_cast<int>(control_points.size()) - 1;
            std::vector<vec4> CK(1);

            curve_derivs_alg1(n, _degree, knots, control_points, t, 1, CK, _periodic);
            auto der = CK[1];

            //printf("(%f,%f,%f,%f)\n",der.x,der.y,der.z,der.w);
            result.x = der.x;
            result.y = der.y;
            result.z = der.z;
        }


        double length(const double t0, const double t1, const double tol = CMMCORE_NURBS_INTEGRATION_TOL) const {
            double result = 0.0;
            double error = 0.0;

            integrate([this](const double t)-> double {
                          vec3 der;
                          derivative(t, der);
                          return der.length();
                      },
                      t0, t1,
                      result,
                      error, tol);

            return result;
        }

        double length(const double t, const double tol = CMMCORE_NURBS_INTEGRATION_TOL) const {
            return length(_interval[0], t, tol);
        }

        double length() const {
            return length(_interval[0], _interval[1], CMMCORE_NURBS_INTEGRATION_TOL);
        };

        double lengthAt(const double l, const double t0, const double initial_guess,
                        const double tol = CMMCORE_NURBS_INTEGRATION_TOL) {
            vec3 der;
            auto fun = [&](double t)-> double {
                derivative(t, der);
                return der.length();
            };

            return find_t1_newton(t0, l, fun, initial_guess, tol);
        }

        double lengthAt(const double l, const double tol = CMMCORE_NURBS_INTEGRATION_TOL) {
            return lengthAt(l, _interval[0], (_interval[0] + _interval[1]) * 0.5, tol);
        }

        void derivative(const double t, std::vector<vec3> &result, const size_t d) const {
            int n = static_cast<int>(control_points.size()) - 1;
            result.resize(d);

            std::vector<vec4> CK(result.size());
            curve_derivs_alg1(n, _degree, knots, control_points, t, d, CK, _periodic);
            for (int i = 0; i < d; ++i) {
                auto &der = CK[i];
                auto &target = result[i];
                target.x = der.x;
                target.y = der.y;
                target.z = der.z;
            }
        }

        void derivative(const double t, std::vector<vec3> &result) const {
            derivative(t, result, result.size());
        }

        void derivatives(const double t, vec3 &d0, vec3 &d1, vec3 &d2) const {
            vec3 temp;
            evaluate(t, d0);
            evaluate(t - 1e-5, temp);
            evaluate(t + 1e-5, d2);
            d1.set(d2 - temp);
            d1 /= 2e-5;
            d2 -= (2 * d0);
            d2 += temp;
            d2 /= 1e-10;
        }

        void derivative_fdm(const double t, vec3 &result) const {
            vec3 temp;
            evaluate(t - 1e-5, temp);
            evaluate(t + 1e-5, result);
            result -= temp;
            result /= 2e-5;
        }

        /**
         ..math:
         \NURBSCurve
         \frac{d^2 f(t)}{dt^2} \approx \frac{f(t + h) - 2f(t) + f(t - h)}{h^2}
         */
        void second_derivative_fdm(const double t, vec3 &result) const {
            vec3 temp;
            evaluate(t - 1e-5, temp);
            evaluate(t + 1e-5, result);
            result -= temp;
            result /= 2e-5;
        }


        vec3 &operator()(const double t, vec3 &result) const {
            evaluate(t, result);
            return result;
        }

        vec3 operator()(const double t) const {
            vec3 result(0);
            evaluate(t, result);
            return result;
        }

        void get_derivative(NURBSCurve &derivative) const {
            get_nurbs_derivative(_degree, knots, control_points, derivative._degree, derivative.knots,
                                 derivative.control_points);

            derivative.update_interval();
        }

        NURBSCurve get_derivative() const {
            NURBSCurve derivative{};
            get_nurbs_derivative(_degree, knots, control_points, derivative._degree, derivative.knots,
                                 derivative.control_points);
            derivative.update_interval();
            return derivative;
        }


        const std::vector<vec4> &get_control_points() const {
            return control_points;
        }

        std::vector<vec3> &get_control_points3d(std::vector<vec3> &result) const {
            result.resize(control_points.size());
            for (size_t i = 0; i < control_points.size(); ++i) {
                result[i].set(control_points[i].to_vec3());
            }
            return result;
        }

        std::vector<vec3> get_control_points3d() const {
            std::vector<vec3> result;
            result.resize(control_points.size());
            for (size_t i = 0; i < control_points.size(); ++i) {
                result[i].set(control_points[i].to_vec3());
            }
            return result;
        }

        void set_control_points(const std::vector<vec4> &cpts) {
            bool change_size = (cpts.size() != control_points.size());
            control_points = std::move(cpts);
            if (change_size) {
                if (_periodic) {
                    generate_knots_periodic();
                } else {
                    generate_knots();
                }
            }
        }

        int get_degree() const {
            return _degree;
        }

        void set_degree(const int val) {
            _degree = val;
            update_interval();
        }

        const std::vector<double> &get_knots() const {
            return knots;
        }

        void set_knots(std::vector<double> &knts) {
            knots = std::move(knts);
            update_interval();
        }

        const std::array<double, 2> &interval() const {
            return _interval;
        }

        const AABB &aabb() const {
            return _aabb;
        }

        // Member variables
        std::vector<vec4> control_points{}; // Each control point is (x, y, z, weight)
        int _degree = 0;
        std::vector<double> knots{};
        bool _periodic = false;
        std::vector<vec3> _hull{};
        std::array<double, 2> _interval{0., 1.}; // [min_knot, max_knot]
        AABB _aabb{};
        // Helper methods
        void rebuildAABB() {
            _aabb.min.set(control_points[0].to_vec3());
            _aabb.max.set(control_points[0].to_vec3());
            for (auto p: control_points) {
                _aabb.expand(p.to_vec3());
            }
        }

        AABB &bbox() {
            rebuildAABB();
            return _aabb;
        }

        // Function to decompose a NURBS curve into Bezier curves
        void decompose(std::vector<NURBSCurve> &bezierCurves) const {
            // Ensure that the curve is valid for decomposition
            //assert(!curve.knots.empty());
            //assert(curve._degree >= 1);


            // Create a temporary copy of the curve to perform splits
            NURBSCurve temp = *this;
            temp.update_interval();


            // Extract unique knots from the knot vector
            std::vector<double> unique_knots;
            uniqueKnots(temp.knots, unique_knots);

            // Iterate over the unique knots, excluding the first and last
            for (size_t i = 1; i < unique_knots.size() - 1; ++i) {
                double knot = unique_knots[i];

                // Split the curve at the current knot
                auto [first, second] = temp.split(knot);

                // Add the first segment to the Bezier curves list
                bezierCurves.push_back(first);

                // Continue with the second segment for further splitting
                temp = second;
            }

            // Add the final segment to the Bezier curves list
            bezierCurves.push_back(temp);
        }

        /**
             * Join two NURBS curves if they share an endpoint
             * @param first First curve to join
             * @param second Second curve to join
             * @param result The resulting joined curve
             * @return true if curves were successfully joined, false otherwise
             */
        bool join(const NURBSCurve &second, NURBSCurve &result) {
            // Check if either curve is already closed
            if (is_periodic() || second.is_periodic()) {
                return false;
            }

            // Get endpoints of both curves
            vec3 first_start, first_end, second_start, second_end;
            evaluate(_interval[0], first_start);
            evaluate(_interval[1], first_end);
            second.evaluate(second._interval[0], second_start);
            second.evaluate(second._interval[1], second_end);

            // Tolerance for endpoint matching
            constexpr double tol = 1e-10;

            // Check all possible connection combinations
            bool start_start = first_start.distance(second_start) < tol;
            bool start_end = first_start.distance(second_end) < tol;
            bool end_start = first_end.distance(second_start) < tol;
            bool end_end = first_end.distance(second_end) < tol;

            // If no endpoints match, cannot join
            if (!start_start && !start_end && !end_start && !end_end) {
                return false;
            }

            // Get control points from both curves
            std::vector<vec4> new_control_points;
            new_control_points.reserve(control_points.size() + second.control_points.size());

            // Join based on which endpoints match
            if (end_start) {
                // Normal case: first end to second start
                new_control_points = control_points;
                new_control_points.insert(new_control_points.end(),
                                          second.control_points.begin(), second.control_points.end());
            } else if (start_start) {
                // Reverse first curve
                new_control_points.insert(new_control_points.end(),
                                          control_points.rbegin(), control_points.rend());
                new_control_points.insert(new_control_points.end(),
                                          second.control_points.begin(), second.control_points.end());
            } else if (end_end) {
                // Reverse second curve
                new_control_points = control_points;
                new_control_points.insert(new_control_points.end(),
                                          second.control_points.rbegin(), second.control_points.rend());
            } else if (start_end) {
                // Reverse first curve and reverse second curve
                new_control_points.insert(new_control_points.end(),
                                          control_points.rbegin(), control_points.rend());
                new_control_points.insert(new_control_points.end(),
                                          second.control_points.rbegin(), second.control_points.rend());
            }

            // If both start and end points match (after potential reversals), make it periodic
            bool should_be_periodic = false; {
                vec3 new_start, new_end;
                vec4 start_pt = new_control_points.front();
                vec4 end_pt = new_control_points.back();
                new_start.set(start_pt.x / start_pt.w, start_pt.y / start_pt.w, start_pt.z / start_pt.w);
                new_end.set(end_pt.x / end_pt.w, end_pt.y / end_pt.w, end_pt.z / end_pt.w);
                should_be_periodic = new_start.distance(new_end) < tol;
            }

            // Create result curve
            result = NURBSCurve(new_control_points, _degree, should_be_periodic);
            return true;
        }

             /**
         * Writes the NURBS Curve as a STEP entity.
         * Includes weights and other necessary attributes for AP242 compliance.
         *
         * @param out Output stream to write the STEP entities.
         * @param entity_id Reference to the current entity ID, incremented as entities are written.
         * @return The STEP entity ID of the created B-spline curve.
         */
int toSTEP(std::ostream &out, int &entity_id) const {
            // Write CARTESIAN_POINT entities with homogeneous coordinates
            std::vector<int> cp_ids(control_points.size());
            for (size_t i = 0; i < control_points.size(); ++i) {
                const vec4 &cp = control_points[i];
                cp_ids[i] = entity_id;
                out << "#" << entity_id++ << " = CARTESIAN_POINT('', ("
                    << cp.x << ", " << cp.y << ", " << cp.z << "));\n";
            }

            // Get knot multiplicities and unique knots
            std::vector<double> unique_knots;
            std::vector<int> multiplicities;
            getKnotMultiplicities(knots, unique_knots, multiplicities);

            // Convert weights vector to string
            std::stringstream weights_stream;
            weights_stream << "(";
            for (size_t i = 0; i < control_points.size(); ++i) {
                weights_stream << control_points[i].w;
                if (i < control_points.size() - 1) weights_stream << ", ";
            }
            weights_stream << ")";

            // Write the RATIONAL_B_SPLINE_CURVE_WITH_KNOTS entity
            int curve_id = entity_id++;
            out << "#" << curve_id << " = RATIONAL_B_SPLINE_CURVE_WITH_KNOTS('', "
                    << _degree << ", (";

            // Write control point references
            for (size_t i = 0; i < cp_ids.size(); ++i) {
                out << "#" << cp_ids[i];
                if (i < cp_ids.size() - 1) out << ", ";
            }

            out << "), .UNSPECIFIED., .F., .F., "
                    << vectorToString(multiplicities) << ", "
                    << vectorToString(unique_knots) << ", "
                    << ".UNSPECIFIED., "
                    << weights_stream.str()
                    << ");\n";

            return curve_id;
        }
    };

    // External function declarations


    class NURBSSurface {
    public:
        AABB _aabb{};

        NURBSSurface() = default;

        NURBSSurface(const std::vector<std::vector<vec4> > &control_points, const std::array<int, 2> &degree,
                     const std::vector<double> &knots_u = {}, const std::vector<double> &knots_v = {})
            : _degree(degree), _control_points(control_points) {
            if (control_points.empty() || control_points[0].empty()) {
                throw std::invalid_argument("Control points cannot be empty");
            }
            _size = {control_points.size(), control_points[0].size()};
            if (knots_u.empty()) {
                generate_knots_u();
            } else {
                _knots_u = knots_u;
            }
            if (knots_v.empty()) {
                generate_knots_v();
            } else {
                _knots_v = knots_v;
            }
            update_interval();
        }

        NURBSSurface(const std::vector<std::vector<vec3> > &control_points, const std::array<int, 2> &degree,
                     const std::vector<double> &knots_u = {}, const std::vector<double> &knots_v = {}
        ): _degree(degree) {
            _control_points.resize(control_points.size());
            for (size_t i = 0; i < control_points.size(); ++i) {
                _control_points[i].resize(control_points[i].size());

                for (size_t j = 0; j < control_points[i].size(); ++j) {
                    _control_points[i][j].set(control_points[i][j]);
                    _control_points[i][j].w = 1.;
                }
            }
            _size = {control_points.size(), control_points[0].size()};
            if (knots_u.empty()) {
                generate_knots_u();
            } else {
                _knots_u = knots_u;
            }
            if (knots_v.empty()) {
                generate_knots_v();
            } else {
                _knots_v = knots_v;
            }
            update_interval();
        }

        std::array<size_t, 2> &shape() {
            return _size;
        }

        std::array<size_t, 2> shape() const {
            return _size;
        }

        std::vector<double> &knots_u() {
            return _knots_u;
        }

        std::vector<double> knots_u() const {
            return _knots_u;
        }

        std::vector<double> &knots_v() {
            return _knots_v;
        }

        std::vector<double> knots_v() const {
            return _knots_v;
        }

        std::array<std::array<double, 2>, 2> &interval() {
            return _interval;
        }

        std::array<std::array<double, 2>, 2> interval() const {
            return _interval;
        }

        void generate_knots_u() {
            std::size_t nu = _size[0];
            _knots_u.resize(nu + _degree[0] + 1);
            std::fill(_knots_u.begin(), _knots_u.begin() + _degree[0] + 1, 0.0);
            std::iota(_knots_u.begin() + _degree[0] + 1, _knots_u.end() - _degree[0], 1.0);
            std::fill(_knots_u.end() - _degree[0], _knots_u.end(), static_cast<double>(nu - _degree[0]));
            update_interval_u();
        }

        void generate_knots_v() {
            std::size_t nv = _size[1];
            _knots_v.resize(nv + _degree[1] + 1);
            std::fill(_knots_v.begin(), _knots_v.begin() + _degree[1] + 1, 0.0);
            std::iota(_knots_v.begin() + _degree[1] + 1, _knots_v.end() - _degree[1], 1.0);
            std::fill(_knots_v.end() - _degree[1], _knots_v.end(), static_cast<double>(nv - _degree[1]));
            update_interval_v();
        }

        void evaluate(double u, double v, vec3 &result) const noexcept {
            surface_point(_size[0] - 1, _degree[0], _knots_u, _size[1] - 1, _degree[1], _knots_v, _control_points, u, v,
                          0, 0, result);
        }

        void derivative_u(double u, double v, vec3 &result) const {
            vec3 temp1, temp2;
            evaluate(u - 1e-5, v, temp1);
            evaluate(u + 1e-5, v, temp2);
            temp2 -= temp1;
            result.set(temp2 / (2 * 1e-5
                       ));
        }

        void derivative_v(double u, double v, vec3 &result) const {
            vec3 temp1, temp2;
            evaluate(u, v - 1e-5, temp1);
            evaluate(u, v + 1e-5, temp2);
            temp2 -= temp1;
            result.set(temp2 / (2 * 1e-5
                       ));
        }

        void normal(double u, double v, vec3 &result) const {
            vec3 temp1, temp2;
            derivative_u(u, v, temp1);
            derivative_v(u, v, temp2);
            temp1.cross(temp2, result);
            result.normalize();
        }

        void insert_knot(const double t, const int direction, const int count) {
            if (direction == 0) {
                int s = find_multiplicity(t, _knots_u);
                int span = find_span(_control_points.size() - 1, _degree[0], t, _knots_u, false);

                auto kv = knot_insertion_kv(_knots_u, t, span, count);
                std::vector<std::vector<vec4> > cpts_tmp;
                std::vector<vec4> cpts = control_points_flat();

                for (size_t v = 0; v < _size[1]; ++v) {
                    std::vector<vec4> ccu;
                    for (size_t u = 0; u < _size[0]; ++u) {
                        ccu.push_back(cpts[u * _size[1] + v]);
                    }


                    cpts_tmp.push_back(knot_insertion(_degree[0], _knots_u, ccu, t, count, s, span));
                }


                nurbs_helpers::flipControlPointsU(cpts_tmp, _control_points);

                _size[0] += count;

                _knots_u = kv;
            } else {
                insert_knot_v(t, count);
            }
        }

        void insert_knot_u(double t, int count) {
            insert_knot(t, 0, count);
        }

        void insert_knot_v(double t, int count) {
            std::size_t new_count_u = _control_points.size();
            std::size_t new_count_v = _control_points[0].size() + count;
            auto cpts = _control_points;
            int span = find_span(static_cast<int>(_size[1] - 1), _degree[1], t, _knots_v, false);
            auto k_v = knot_insertion_kv(_knots_v, t, span, count);
            int s_v = find_multiplicity(t, _knots_v);
            std::vector<std::vector<vec4> > new_control_points(new_count_u,
                                                               std::vector<vec4>(new_count_v, {0, 0, 0, 1}));
            for (std::size_t u = 0; u < new_count_u; ++u) {
                std::vector<vec4> row(new_count_v);
                knot_insertion(_degree[1], _knots_v, cpts[u], t, count, s_v, span, row);
                new_control_points[u] = row;
            }
            _control_points = std::move(new_control_points);
            _knots_v = std::move(k_v);
            _size[1] = new_count_v;
            update_interval();
        }

        std::pair<NURBSSurface, NURBSSurface> split_surface_u(double param, double tol = 1e-7) const {
            if (param <= _interval[0][0] || param >= _interval[0][1] || std::fabs(param - _interval[0][0]) <= tol ||
                std::fabs(param - _interval[0][1]) <= tol) {
                throw std::invalid_argument("Cannot split from the domain edge");
            }
            int ks = find_span(static_cast<int>(_size[0] - 1), _degree[0], param, _knots_u, false) - _degree[0] + 1;
            int s = find_multiplicity(param, _knots_u);

            int r = _degree[0] - s;
            NURBSSurface temp_obj = *this;

            temp_obj.insert_knot_u(param, r);
            std::vector<double> temp_knots_u = temp_obj._knots_u;
            std::vector<std::vector<vec4> > tcpts = temp_obj._control_points;

            int knot_span = find_span(static_cast<int>(temp_obj._size[0] - 1), _degree[0], param, temp_knots_u,
                                      false) + 1;
            std::vector<double> surf1_kv(temp_obj._knots_u.begin(), temp_obj._knots_u.begin() + knot_span);
            std::vector<double> surf2_kv(temp_obj._knots_u.begin() + knot_span, temp_obj._knots_u.end());
            surf1_kv.push_back(param);
            for (int j = 0; j <= _degree[0]; ++j) {
                surf2_kv.insert(surf2_kv.begin(), param);
            }
            std::vector<std::vector<vec4> > surf1_ctrlpts(temp_obj._control_points.begin(),
                                                          temp_obj._control_points.begin() + ks + r);
            std::vector<std::vector<vec4> > surf2_ctrlpts(temp_obj._control_points.begin() + ks + r - 1,
                                                          temp_obj._control_points.end());
            NURBSSurface surf1(surf1_ctrlpts, {_degree[0], _degree[1]}, surf1_kv, _knots_v);
            NURBSSurface surf2(surf2_ctrlpts, {_degree[0], _degree[1]}, surf2_kv, _knots_v);
            return {surf1, surf2};
        }

        std::pair<NURBSSurface, NURBSSurface> split_surface_v(double param, double tol = 1e-7) const {
            if (param <= _interval[1][0] || param >= _interval[1][1] || std::fabs(param - _interval[1][0]) <= tol ||
                std::fabs(param - _interval[1][1]) <= tol) {
                throw std::invalid_argument("Cannot split from the domain edge");
            }
            int ks = find_span(static_cast<int>(_size[1] - 1), _degree[1], param, _knots_v, false) - _degree[1] + 1;
            int s = find_multiplicity(param, _knots_v);
            int r = _degree[1] - s;
            NURBSSurface temp_obj = *this;
            temp_obj.insert_knot_v(param, r);
            int knot_span = find_span(static_cast<int>(temp_obj._size[1] - 1), _degree[1], param, temp_obj._knots_v,
                                      false) + 1;
            std::vector<double> surf1_kv(temp_obj._knots_v.begin(), temp_obj._knots_v.begin() + knot_span);
            std::vector<double> surf2_kv(temp_obj._knots_v.begin() + knot_span, temp_obj._knots_v.end());
            surf1_kv.push_back(param);
            for (int j = 0; j <= _degree[1]; ++j) {
                surf2_kv.insert(surf2_kv.begin(), param);
            }
            std::vector<std::vector<vec4> > surf1_ctrlpts(_size[0]);
            std::vector<std::vector<vec4> > surf2_ctrlpts(_size[0]);
            for (std::size_t i = 0; i < _size[0]; ++i) {
                surf1_ctrlpts[i].assign(temp_obj._control_points[i].begin(),
                                        temp_obj._control_points[i].begin() + ks + r);
                surf2_ctrlpts[i].assign(temp_obj._control_points[i].begin() + ks + r - 1,
                                        temp_obj._control_points[i].end());
            }
            NURBSSurface surf1(surf1_ctrlpts, {_degree[0], _degree[1]}, _knots_u, surf1_kv);
            NURBSSurface surf2(surf2_ctrlpts, {_degree[0], _degree[1]}, _knots_u, surf2_kv);
            return {surf1, surf2};
        }

        std::vector<vec4> control_points_flat() const {
            std::vector<vec4> control_points_flt(_size[0] * _size[1]);
            for (std::size_t i = 0; i < _size[0]; ++i) {
                for (std::size_t j = 0; j < _size[1]; ++j) {
                    control_points_flt[i * _size[1] + j].set(_control_points[i][j]);
                }
            };
            return control_points_flt;
        }

        void control_points_flat3d(std::vector<vec3> &control_points_flt) const {
            control_points_flt.resize(_size[0] * _size[1]);
            for (std::size_t i = 0; i < _size[0]; ++i) {
                for (std::size_t j = 0; j < _size[1]; ++j) {
                    const double w = _control_points[i][j].w;
                    control_points_flt[i * _size[1] + j].set(_control_points[i][j].x / w, _control_points[i][j].y / w,
                                                             _control_points[i][j].z / w);
                }
            };
        }

        std::vector<vec3> control_points_flat3d() const {
            std::vector<vec3> control_points_flt;
            control_points_flat3d(control_points_flt);
            return control_points_flt;
        }


        const std::vector<std::vector<vec4> > &control_points() const {
            return _control_points;
        }

        std::vector<std::vector<vec3> > control_points3d() const {
            double w;
            std::vector<std::vector<vec3> > control_points_(_size[0], std::vector<vec3>(_size[1]));
            for (std::size_t i = 0; i < _size[0]; ++i) {
                for (std::size_t j = 0; j < _size[1]; ++j) {
                    w = _control_points[i][j].w;
                    control_points_[i][j].set(_control_points[i][j].x / w, _control_points[i][j].y / w,
                                              _control_points[i][j].z / w);
                }
            }
            return control_points_;
        }

        std::vector<std::vector<vec4> > &control_points() {
            return _control_points;
        }

        std::array<NURBSSurface, 4> subdivide(

        ) const {
            double umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            double vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            std::array<NURBSSurface, 4> surfs;
            //auto umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            //auto vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            auto [s1,s2] = split_surface_u(umid);
            auto [_s11,_s12] = s1.split_surface_v(vmid);
            auto [_s21,_s22] = s2.split_surface_v(vmid);
            surfs[0] = std::move(_s11);
            surfs[1] = std::move(_s12);
            surfs[2] = std::move(_s21);
            surfs[3] = std::move(_s22);
            return surfs;
        }

        std::array<NURBSSurface, 4> subdivide(double umid, double vmid

        ) const {
            std::array<NURBSSurface, 4> surfs;
            //auto umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            //auto vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            auto [s1,s2] = split_surface_u(umid);
            auto [_s11,_s12] = s1.split_surface_v(vmid);
            auto [_s21,_s22] = s2.split_surface_v(vmid);
            surfs[0] = std::move(_s11);
            surfs[1] = std::move(_s12);
            surfs[2] = std::move(_s21);
            surfs[3] = std::move(_s22);
            return surfs;
        }

        void subdivide(
            NURBSSurface &s11,
            NURBSSurface &s12,
            NURBSSurface &s21,
            NURBSSurface &s22,
            double umid,
            double vmid
        ) const {
            //auto umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            //auto vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            auto [s1,s2] = split_surface_u(umid);
            auto [_s11,_s12] = s1.split_surface_v(vmid);
            auto [_s21,_s22] = s2.split_surface_v(vmid);
            s11 = std::move(_s11);
            s12 = std::move(_s12);
            s21 = std::move(_s21);
            s22 = std::move(_s22);
        }

        void subdivide(
            NURBSSurface &s11,
            NURBSSurface &s12,
            NURBSSurface &s21,
            NURBSSurface &s22

        ) const {
            auto umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            auto vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            subdivide(s11, s12, s21, s22, umid, vmid);
        }

        const AABB &aabb() const {
            return _aabb;
        }

        AABB &bbox() {
            _bbox.min.set(_control_points[0][0].to_vec3());
            _bbox.max.set(_control_points[0][0].to_vec3());


            for (std::size_t i = 0; i < _size[0]; ++i) {
                for (std::size_t j = 0; j < _size[1]; ++j) {
                    auto &pt = _control_points[i][j];
                    _bbox.expand(pt.to_vec3());
                }
            }

            return _bbox;
        }

        std::array<int, 2> _degree{};
        std::array<std::size_t, 2> _size{};
        std::array<std::array<double, 2>, 2> _interval{};
        std::vector<std::vector<vec4> > _control_points{};
        std::vector<double> _knots_u{};
        std::vector<double> _knots_v{};
        AABB _bbox{
            {0., 0., 0.}, {0., 0., 0.}
        };


        void update_interval_u() {
            _interval[0][0] = *(_knots_u.begin());
            _interval[0][1] = *(_knots_u.end() - 1);
        }

        void update_interval_v() {
            _interval[1][0] = *(_knots_v.begin());
            _interval[1][1] = *(_knots_v.end() - 1);
        }

        void update_interval() {
            _interval[0][0] = *(_knots_u.begin());
            _interval[0][1] = *(_knots_u.end() - 1);
            _interval[1][0] = *(_knots_v.begin());
            _interval[1][1] = *(_knots_v.end() - 1);
        }

        void get_derivative(int direction, NURBSSurface &result) const {
            if (direction != 0 && direction != 1) {
                throw std::invalid_argument("Direction must be 0 (u) or 1 (v)");
            }

            std::vector<std::vector<vec4> > new_control_points;
            std::vector<double> new_knots;
            int new_degree;

            if (direction == 0) {
                // u-direction
                new_degree = _degree[0] - 1;
                new_knots = std::vector<double>(_knots_u.begin() + 1, _knots_u.end() - 1);
                new_control_points.resize(_size[0] - 1, std::vector<vec4>(_size[1]));

                for (std::size_t j = 0; j < _size[1]; ++j) {
                    for (std::size_t i = 0; i < _size[0] - 1; ++i) {
                        double factor = _degree[0] / (_knots_u[i + _degree[0] + 1] - _knots_u[i + 1]);
                        new_control_points[i][j] = (_control_points[i + 1][j] - _control_points[i][j]) * factor;
                    }
                }
            } else {
                // v-direction
                new_degree = _degree[1] - 1;
                new_knots = std::vector<double>(_knots_v.begin() + 1, _knots_v.end() - 1);
                new_control_points.resize(_size[0], std::vector<vec4>(_size[1] - 1));

                for (std::size_t i = 0; i < _size[0]; ++i) {
                    for (std::size_t j = 0; j < _size[1] - 1; ++j) {
                        double factor = _degree[1] / (_knots_v[j + _degree[1] + 1] - _knots_v[j + 1]);
                        new_control_points[i][j] = (_control_points[i][j + 1] - _control_points[i][j]) * factor;
                    }
                }
            }

            std::array<int, 2> new_degrees = {
                direction == 0 ? new_degree : _degree[0], direction == 1 ? new_degree : _degree[1]
            };
            result._control_points = new_control_points;
            result._degree = new_degrees;
            result._knots_u = direction == 0 ? new_knots : _knots_u;
            result._knots_v = direction == 1 ? new_knots : _knots_v;
        }

        NURBSSurface get_derivative(int direction) const {
            NURBSSurface result;
            get_derivative(direction, result);
            return result;
        }


        /**
        * Extract an isocurve from a NURBS surface at a given parameter in the u or v direction
        */
        /**
         * Extracts an isocurve from a NURBS surface at a given parameter in the u or v direction.
         *
         * @param param The parameter value where the isocurve is extracted.
         * @param direction The direction (U or V) in which to extract the isocurve.
         * @return A NURBSCurve object representing the extracted isocurve.
         */
        NURBSCurve extract_isocurve(const double param, const SurfaceParameter direction = SurfaceParameter::U) const {
            if (direction != SurfaceParameter::U && direction != SurfaceParameter::V) {
                throw std::invalid_argument("Direction must be either 'U' or 'V'.");
            }


            if (direction == SurfaceParameter::U) {
                // For u-direction: fix u and vary v
                auto param_range = _interval[0]; // u range
                if (param < param_range[0] || param > param_range[1]) {
                    throw std::invalid_argument("Parameter " + std::to_string(param) +
                                                " is out of range [" + std::to_string(param_range[0]) +
                                                "," + std::to_string(param_range[1]) + "]");
                }

                // Find span and basis functions in u direction
                int n_u = _size[0] - 1;
                int degree_u = _degree[0];
                int span = find_span(n_u, degree_u, param, _knots_u, false);
                std::array<double, CMMCORE_DEG_STACK_SIZE> basis{};
                basis_funs(span, param, degree_u, _knots_u, basis);


                // The curve will have as many control points as the surface has in v direction
                size_t m = _size[1];
                NURBSCurve crv(m, _degree[1]);
                crv.knots = _knots_v;
                crv.update_interval();


                // Compute control points for the extracted curve
                for (int i = 0; i < m; ++i) {
                    vec4 point(0.0, 0.0, 0.0, 0.0);
                    for (int j = 0; j <= degree_u; ++j) {
                        point += _control_points[span - degree_u + j][i] * basis[j];
                    }
                    crv.control_points[i] = point;
                }
                crv.update_interval();

                return crv;
            } else {
                // For v-direction: fix v and vary u
                auto param_range = _interval[1]; // v range
                if (param < param_range[0] || param > param_range[1]) {
                    throw std::invalid_argument("Parameter " + std::to_string(param) +
                                                " is out of range [" + std::to_string(param_range[0]) +
                                                "," + std::to_string(param_range[1]) + "]");
                }

                // Find span and basis functions in v direction
                int n_v = _size[1] - 1;
                int degree_v = _degree[1];
                int span = find_span(n_v, degree_v, param, _knots_v, false);
                std::array<double, CMMCORE_DEG_STACK_SIZE> basis{};
                basis_funs(span, param, degree_v, _knots_v, basis);


                // The curve will have as many control points as the surface has in u direction
                int m = _size[0];
                NURBSCurve crv(m, _degree[0]);
                crv.knots = _knots_u;
                crv.update_interval();

                // Compute control points for the extracted curve
                for (int i = 0; i < m; ++i) {
                    vec4 point(0.0, 0.0, 0.0, 0.0);
                    for (int j = 0; j <= degree_v; ++j) {
                        point += _control_points[i][span - degree_v + j] * basis[j];
                    }
                    crv.control_points[i] = point;
                }

                return crv;
            }
        }


        /**
         * Extracts the boundary NURBS curves from the surface at the given parameters.
         *
         * @param result Output array of 4 NURBSCurve objects representing the boundaries.
         */
        void extract_surface_boundaries(std::array<NURBSCurve, 4> &result) const {
            const double u_min = _interval[0][0];
            const double u_max = _interval[0][1];
            const double v_min = _interval[1][0];
            const double v_max = _interval[1][1];

            // Extract iso-curves at the boundaries
            result[0] = extract_isocurve(u_min, SurfaceParameter::U); // v-direction curve at u_min
            result[1] = extract_isocurve(u_max, SurfaceParameter::U); // v-direction curve at u_max
            result[2] = extract_isocurve(v_min, SurfaceParameter::V); // u-direction curve at v_min
            result[3] = extract_isocurve(v_max, SurfaceParameter::V); // u-direction curve at v_max
        }

       /**
         * Writes comprehensive STEP entities for the NURBS surface, including geometry and presentation.
         * Updated to comply with AP242 schema.
         */
        int toSTEP(std::ostream &out, int &entity_id) const {
    // Write CARTESIAN_POINT entities with normalized coordinates
    std::vector<std::vector<int>> cp_ids(_size[0], std::vector<int>(_size[1]));
    for (std::size_t i = 0; i < _size[0]; ++i) {
        for (std::size_t j = 0; j < _size[1]; ++j) {
            const vec4 &cp = _control_points[i][j];
            // Normalize coordinates by weight
            double x = cp.x / cp.w;
            double y = cp.y / cp.w;
            double z = cp.z / cp.w;
            int id = entity_id++;
            cp_ids[i][j] = id;
            out << "#" << id << " = CARTESIAN_POINT('', (" << x << ", " << y << ", " << z << "));\n";
        }
    }

    // Prepare control points grid string
    std::stringstream cp_grid_stream;
    cp_grid_stream << "(";
    for (std::size_t i = 0; i < _size[0]; ++i) {
        cp_grid_stream << "(";
        for (std::size_t j = 0; j < _size[1]; ++j) {
            cp_grid_stream << "#" << cp_ids[i][j];
            if (j < _size[1] - 1)
                cp_grid_stream << ", ";
        }
        cp_grid_stream << ")";
        if (i < _size[0] - 1)
            cp_grid_stream << ", ";
    }
    cp_grid_stream << ")";

    // Prepare unique knot vectors and multiplicities
    std::vector<double> u_unique_knots, v_unique_knots;
    std::vector<int> u_multiplicities, v_multiplicities;
    getKnotMultiplicities(_knots_u, u_unique_knots, u_multiplicities);
    getKnotMultiplicities(_knots_v, v_unique_knots, v_multiplicities);

    // Convert vectors to STEP-compatible strings
    std::string u_knots_str = vectorToString(u_unique_knots);
    std::string v_knots_str = vectorToString(v_unique_knots);
    std::string u_mult_str = vectorToString(u_multiplicities);
    std::string v_mult_str = vectorToString(v_multiplicities);

    // Prepare weights grid string
    std::stringstream weights_stream;
    weights_stream << "(";
    for (std::size_t i = 0; i < _size[0]; ++i) {
        weights_stream << "(";
        for (std::size_t j = 0; j < _size[1]; ++j) {
            double w = _control_points[i][j].w;
            weights_stream << w;
            if (j < _size[1] - 1)
                weights_stream << ", ";
        }
        weights_stream << ")";
        if (i < _size[0] - 1)
            weights_stream << ", ";
    }
    weights_stream << ")";

    // Write RATIONAL_B_SPLINE_SURFACE_WITH_KNOTS entity
    int bspline_surface_id = entity_id++;
    out << "#" << bspline_surface_id << " = RATIONAL_B_SPLINE_SURFACE_WITH_KNOTS('', "
        << _degree[0] << ", " << _degree[1] << ", "
        << cp_grid_stream.str() << ", "
        << ".UNSPECIFIED., .F., .F., .F., "
        << u_mult_str << ", " << u_knots_str << ", "
        << v_mult_str << ", " << v_knots_str << ", "
        << ".UNSPECIFIED., "
        << weights_stream.str() << ");\n";
    std::array<NURBSCurve, 4> boundaries;
            extract_surface_boundaries(boundaries) ;

    // Write boundary curves
    std::array<int, 4> curve_ids;
    for (int i = 0; i < 4; ++i) {
        curve_ids[i] = boundaries[i].toSTEP(out, entity_id);
    }

    // Create and write VERTEX_POINT entities for corners
    std::array<int, 4> vertex_ids;
    for (int i = 0; i < 4; ++i) {
        vec3 point;
        double param = (i < 2) ? boundaries[i].interval()[0] : boundaries[i].interval()[1];
        boundaries[i].evaluate(param, point);
        int cartesian_id = entity_id++;
        out << "#" << cartesian_id << " = CARTESIAN_POINT('', (" << point.x << ", " << point.y << ", " << point.z << "));\n";
        vertex_ids[i] = entity_id++;
        out << "#" << vertex_ids[i] << " = VERTEX_POINT('', #" << cartesian_id << ");\n";
    }

    // Create EDGE_CURVE entities
    std::array<int, 4> edge_ids;
    for (int i = 0; i < 4; ++i) {
        int next = (i + 1) % 4;
        int edge_id = entity_id++;
        out << "#" << edge_id << " = EDGE_CURVE('', #" << vertex_ids[i] << ", #" << vertex_ids[next]
            << ", #" << curve_ids[i] << ", .T.);\n";
        edge_ids[i] = edge_id;
    }

    // Create EDGE_LOOP entity
    int edge_loop_id = entity_id++;
    out << "#" << edge_loop_id << " = EDGE_LOOP('', (";
    for (int i = 0; i < 4; ++i) {
        out << "#" << edge_ids[i];
        if (i < 3) out << ", ";
    }
    out << "));\n";

    // Create FACE_OUTER_BOUND entity
    int face_bound_id = entity_id++;
    out << "#" << face_bound_id << " = FACE_OUTER_BOUND('', #" << edge_loop_id << ", .T.);\n";

    // Create ADVANCED_FACE entity
    int face_id = entity_id++;
    out << "#" << face_id << " = ADVANCED_FACE('', (#" << face_bound_id << "), #" << bspline_surface_id << ", .T.);\n";

    // Create OPEN_SHELL entity
    int shell_id = entity_id++;
    out << "#" << shell_id << " = OPEN_SHELL('', (#" << face_id << "));\n";

    // Create SHELL_BASED_SURFACE_MODEL entity
    int shell_based_model_id = entity_id++;
    out << "#" << shell_based_model_id << " = SHELL_BASED_SURFACE_MODEL('shell_1', (#" << shell_id << "));\n";

    // Create styling entities
    int color_rgb_id = entity_id++;
    out << "#" << color_rgb_id << " = COLOUR_RGB('', 0., 0., 0.);\n";

    int fill_area_style_colour_id = entity_id++;
    out << "#" << fill_area_style_colour_id << " = FILL_AREA_STYLE_COLOUR('', #" << color_rgb_id << ");\n";

    int fill_area_style_id = entity_id++;
    out << "#" << fill_area_style_id << " = FILL_AREA_STYLE('', (#" << fill_area_style_colour_id << "));\n";

    int surface_style_fill_area_id = entity_id++;
    out << "#" << surface_style_fill_area_id << " = SURFACE_STYLE_FILL_AREA(#" << fill_area_style_id << ");\n";

    int surface_side_style_id = entity_id++;
    out << "#" << surface_side_style_id << " = SURFACE_SIDE_STYLE('', (#" << surface_style_fill_area_id << "), .F.);\n";

    int surface_style_usage_id = entity_id++;
    out << "#" << surface_style_usage_id << " = SURFACE_STYLE_USAGE(.BOTH., #" << surface_side_style_id << ");\n";

    // Create MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION entity
    int geom_representation_id = entity_id++;
    out << "#" << geom_representation_id << " = MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('', (#" << face_id << "), #" << bspline_surface_id << ");\n";

    // Create STYLED_ITEM entity
    int styled_item_id = entity_id++;
    out << "#" << styled_item_id << " = STYLED_ITEM('', (#" << surface_style_usage_id << "), #" << geom_representation_id << ");\n";

    // Create PRESENTATION_STYLE_ASSIGNMENT entity
    int presentation_style_assignment_id = entity_id++;
    out << "#" << presentation_style_assignment_id << " = PRESENTATION_STYLE_ASSIGNMENT((#" << styled_item_id << "));\n";

    // Create PRESENTATION_LAYER_ASSIGNMENT entity
    int presentation_layer_assignment_id = entity_id++;
    out << "#" << presentation_layer_assignment_id << " = PRESENTATION_LAYER_ASSIGNMENT('Default','', (#" << styled_item_id << "));\n";

    // Create PRODUCT entity
    int product_id = entity_id++;
    out << "#" << product_id << " = PRODUCT('Document','Document','', (#" << presentation_layer_assignment_id << "));\n";

    // Create PRODUCT_DEFINITION_CONTEXT entity
    int product_definition_context_id = entity_id++;
    out << "#" << product_definition_context_id << " = PRODUCT_DEFINITION_CONTEXT('3D Mechanical Parts', #" << (entity_id + 1) << ", 'design');\n";

    // Create PRODUCT_DEFINITION entity
    int product_definition_id = entity_id++;
    out << "#" << product_definition_id << " = PRODUCT_DEFINITION('A','First version', #" << (entity_id + 1) << ", #" << product_definition_context_id << ");\n";

    // Create PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE entity
    int product_definition_formation_id = entity_id++;
    out << "#" << product_definition_formation_id << " = PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('A', 'First version', #" << (entity_id + 1) << ", .MADE.);\n";

    // Create PRODUCT_RELATED_PRODUCT_CATEGORY entity
    int product_related_category_id = entity_id++;
    out << "#" << product_related_category_id << " = PRODUCT_RELATED_PRODUCT_CATEGORY('tool','tool', (#" << product_definition_formation_id << ");\n";

    // Create APPLICATION_PROTOCOL_DEFINITION entity
    int application_protocol_def_id = entity_id++;
    out << "#" << application_protocol_def_id << " = APPLICATION_PROTOCOL_DEFINITION('Draft International Standard', 'automotive_design', 1999, #" << (entity_id + 1) << ");\n";

    // Create APPLICATION_CONTEXT entity
    int application_context_id = entity_id++;
    out << "#" << application_context_id << " = APPLICATION_CONTEXT('data for automotive mechanical design processes');\n";

    // Create PRODUCT_CONTEXT entity
    int product_context_id = entity_id++;
    out << "#" << product_context_id << " = PRODUCT_CONTEXT('3D Mechanical Parts', #" << application_context_id << ", 'mechanical');\n";

    // Create final PRODUCT entity
    int product_id_final = entity_id++;
    out << "#" << product_id_final << " = PRODUCT('Document','Document','', (#" << product_context_id << "));\n";

    // Create SHAPE_REPRESENTATION_RELATIONSHIP entity
    int shape_rep_relationship_id = entity_id++;
    out << "#" << shape_rep_relationship_id << " = SHAPE_REPRESENTATION_RELATIONSHIP('', '', #" << geom_representation_id << ", #" << product_id_final << ");\n";

    return bspline_surface_id;
}

        /**
         * Writes a complete STEP file containing the NURBS surface with comprehensive representations.
         * Modified to comply with AP242 schema.
         */
        void writeSTEPFile(const std::string &filename) const {
            std::ofstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open file for writing");
            }

            // Write the STEP file header
            file << "ISO-10303-21;\n";
            file << "HEADER;\n";
            file << "/* Generated by software containing ST-Developer\n"
                    " * from STEP Tools, Inc. (www.steptools.com)\n"
                    " */\n";
            file << "/* OPTION: using custom schema-name function */\n\n";
            file << "FILE_DESCRIPTION((''), '2;1');\n\n";
            file << "FILE_NAME(\n"
                    "/* name */ 'nc',\n"
                    "/* time_stamp */ '" << getCurrentDateTime()
                 << "',\n"
                    "/* author */ ('Author'),\n"
                    "/* organization */ ('Company'),\n"
                    "/* preprocessor_version */ 'ST-DEVELOPER v19.2',\n"
                    "/* originating_system */ 'Rhino 8.10',\n"
                    "/* authorisation */ 'Authorisation');\n\n";
            file << "FILE_SCHEMA (('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 3 1 4 }'));\n";
            file << "ENDSEC;\n";
            file << "DATA;\n";

            // Initialize entity ID starting from 10
            int entity_id = 10;

            // Write the STEP entities
            toSTEP(file, entity_id);

            // Close the DATA section and file
            file << "ENDSEC;\n";
            file << "END-ISO-10303-21;\n";

            file.close();
        }

        /**
         * Retrieves the current date and time in STEP-compatible format.
         * Example: '2024-11-15T02:15:59+03:00'
         */
        static std::string getCurrentDateTime() {
            std::time_t t = std::time(nullptr);
            std::tm tm = *std::gmtime(&t);
            char buffer[30];
            std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", &tm);
            return std::string(buffer);
        }


        void decompose(const SurfaceParameter direction, std::vector<NURBSSurface> &bezierSurfaces) const {
            //std::vector<double>::iterator knots_it =iterate_unique_knots(_knots_u,_degree[0]);

            NURBSSurface temp = *this;
            std::vector<double> uknots;
            temp.update_interval();
            if (direction == SurfaceParameter::U) {
                uniqueKnots(temp._knots_u, uknots);

                for (auto &knot: uknots) {
                    if (knot != uknots[0] && knot != (uknots[uknots.size() - 1])) {
                        auto [first, second] = temp.split_surface_u(knot);

                        bezierSurfaces.push_back(first);
                        temp = second;
                    }
                }
                bezierSurfaces.push_back(temp);
            } else if (direction == SurfaceParameter::V) {
                uniqueKnots(temp._knots_v, uknots);

                for (auto &knot: uknots) {
                    if (knot != uknots[0] && knot != (uknots[uknots.size() - 1])) {
                        auto [first, second] = temp.split_surface_v(knot);

                        bezierSurfaces.push_back(first);
                        temp = second;
                    }
                }
                bezierSurfaces.push_back(temp);
            }
        }

        void decompose(std::vector<NURBSSurface> &bezierSurfaces) const {
            std::vector<NURBSSurface> bezierSurfaces1;
            decompose(SurfaceParameter::U, bezierSurfaces1);
            for (auto &item: bezierSurfaces1) {
                item.decompose(SurfaceParameter::V, bezierSurfaces);
            }
        }
    };

    inline void surfaceCrossProduct(const NURBSSurface &a, const NURBSSurface &b, NURBSSurface &result) {
        //if (a._size != b._size || a._knots_u != b._knots_u || a._knots_v != b._knots_v) {
        //    throw std::invalid_argument("Surfaces must have the same size and knot vectors");
        //}
        result._control_points.resize(a._size[0],
                                      std::vector<vec4>(a._size[1],
                                                        {0., 0., 0., 1.}
                                      )
        );
        result._knots_u = a._knots_u;
        result._knots_v = b._knots_v;
        result._degree = {a._degree[0], b._degree[1]};
        for (int k = 0; k < a._size[0]; ++k) {
            for (int j = 0; j < a._size[1]; ++j) {
                a._control_points[k][j].cross(b._control_points[j][k], result._control_points[k][j]);
            }
        }
        result.update_interval();
    }

    inline void cross(const NURBSSurface &a, const NURBSSurface &b, NURBSSurface &result) {
        surfaceCrossProduct(a, b, result);
    }

    inline NURBSSurface cross(const NURBSSurface &a, const NURBSSurface &b) {
        NURBSSurface result;
        surfaceCrossProduct(a, b, result);
        return result;
    }
}


#endif //NURBS_H
