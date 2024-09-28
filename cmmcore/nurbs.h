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
#include <iterator>
#include <ranges>
#include <cassert>
#include "cmmcore/nurbs_utils.h"
#include "cmmcore/bvh.hpp"


namespace cmmcore {
    class NURBSCurve {
    public:
        NURBSCurve() = default;

        NURBSCurve(const std::vector<vec4> &control_points): NURBSCurve(
            control_points, control_points.size() >= 4 ? 3 : 1, false) {
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
        NURBSCurve(const std::vector<vec3> &control_points, const int degree,
                   const bool periodic): _degree(degree), _periodic(periodic) {
            // Initialize control points
            const int nnn = control_points.size();
            _control_points.resize(nnn);
            for (size_t i = 0; i < control_points.size(); ++i) {
                _control_points[i].set(control_points[i]);
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
                   bool periodic) : _control_points(control_points), _degree(degree), _periodic(periodic) {
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
            : _control_points(control_points), _degree(degree), _knots(knots), _periodic(periodic) {
            _update_interval();
            rebuildAABB();
        }

        // Copy constructor and assignment operator
        NURBSCurve(const NURBSCurve &other) {
            _control_points = other._control_points;
            _degree = other._degree;
            //for (int i = 0; i < other._knots.size(); ++i) {
            //    printf("%f ",other._knots[i]);
            //}
            //printf("\n ");
            _knots = other._knots;
            //for (int i = 0; i < _knots.size(); ++i) {
            //    printf("%f ",_knots[i]);
            //}
            _aabb = other._aabb;
            _periodic = other._periodic;
            _interval = other._interval;
            _hull = other._hull;
        };

        NURBSCurve &operator=(const NURBSCurve &other) {
            _control_points = other._control_points;
            _degree = other._degree;
            _knots = other._knots;
            _aabb = other._aabb;
            _periodic = other._periodic;
            _interval = other._interval;
            _hull = other._hull;
            return *this;
        };

        // Methods

        // Check if the NURBS curve is periodic
        bool is_periodic() const {
            for (int j = 0; j < _degree; ++j) {
                for (int i = 0; i < 4; ++i) {
                    if (_control_points[j][i] != _control_points[_control_points.size() - _degree + j][i]) {
                        return false;
                    }
                }
            }
            return true;
        }

        // Update interval based on knots
        void _update_interval() {
            _interval[0] = std::min_element(_knots.begin(), _knots.end())[0];
            _interval[1] = std::max_element(_knots.begin(), _knots.end())[0];
        }

        // Hook function called after knots are updated
        void knots_update_hook() {
            _update_interval();
        }

        // Normalize the knot vector
        void normalize_knots() {
            double start = std::min_element(_knots.begin(), _knots.end())[0];
            double end = std::max_element(_knots.begin(), _knots.end())[0];
            double d = 1.0 / (end - start);

            for (auto &knot: _knots) {
                knot = (knot - start) * d;
            }
            knots_update_hook();
        }

        // Generate default open knot vector
        void generate_knots() {
            int n = static_cast<int>(_control_points.size());
            _knots.clear();
            _knots.resize(n + _degree + 1);

            std::fill(_knots.begin(), _knots.begin() + _degree + 1, 0.0);
            std::iota(_knots.begin() + _degree + 1, _knots.end() - _degree, 1.0);
            std::fill(_knots.end() - _degree, _knots.end(), static_cast<double>(n - _degree));

            knots_update_hook();
        }

        // Make the NURBS curve periodic
        void make_periodic() {
            if (is_periodic()) {
                return;
            }
            int n = static_cast<int>(_control_points.size());
            int new_n = n + _degree;
            std::vector<vec4> new_control_points(new_n);

            // Copy original control points
            std::copy(_control_points.begin(), _control_points.end(), new_control_points.begin());

            // Add the first degree control points to the end
            for (int i = 0; i < _degree; ++i) {
                new_control_points[n + i] = _control_points[i];
            }

            _control_points = std::move(new_control_points);
            generate_knots_periodic();
            _periodic = true;
        }

        // Generate knot vector for a periodic NURBS curve
        void generate_knots_periodic() {
            int n = static_cast<int>(_control_points.size());
            int m = n + _degree + 1;
            _knots.resize(m);
            for (int i = 0; i < m; ++i) {
                _knots[i] = static_cast<double>(i - _degree);
            }
            knots_update_hook();
        }

        // Insert knots into the curve
        inline void insert_knot(double t, int count) {
            int n = static_cast<int>(_control_points.size());
            int new_count = n + count;
            const std::vector<vec4> cpts = _control_points;

            // Find knot span
            int span = find_span(n - 1, _degree, t, _knots, false);


            // Compute new knot vector
            std::vector<double> k_v = knot_insertion_kv(_knots, t, span, count);
            int s_u = find_multiplicity(t, _knots);
            // Compute new control points
            _control_points.resize(new_count);

            for (int i = 0; i <= count; ++i) {
                _control_points[i] = cpts.back() - i;
            }
            knot_insertion(_degree, _knots, cpts, t, count, s_u, span, _control_points);


            // Update curve
            _knots = std::move(k_v);
            _update_interval();
        }

        // Split the curve at a given parameter value
        inline std::pair<NURBSCurve, NURBSCurve> split(const double param, const bool normalize_knots = false) const {
            if (param <= _interval[0] || param >= _interval[1] ||
                std::fabs(param - _interval[0]) <= 1e-12 || std::fabs(param - _interval[1]) <= 1e-12) {
                throw std::invalid_argument("Cannot split from the domain edge.");
            }

            int n_ctrlpts = static_cast<int>(_control_points.size());
            auto ks = find_span(n_ctrlpts, _degree, param, _knots, false) - _degree + 1;
            int s = find_multiplicity(param, _knots);
            int r = _degree - s;

            // Insert knot
            NURBSCurve temp_obj(this->_control_points, this->_degree, this->_knots, this->_periodic);
            temp_obj._update_interval();
            temp_obj.insert_knot(param, r);


            // Knot span index
            int knot_span = find_span(temp_obj._control_points.size(), _degree, param, temp_obj._knots, false) + 1;

            // Create knot vectors for the two new curves
            std::vector<double> surf1_kv(knot_span);
            // printf("%d\n",knot_span);
            std::vector<double> surf2_kv(temp_obj._knots.size() - knot_span);

            //surf1_kv.assign(temp_knots.begin(), temp_knots.begin()+ knot_span);
            //surf2_kv.assign(temp_knots.end()-knot_span, temp_knots.end());

            for (int i = 0; i < knot_span; ++i) {
                surf1_kv[i] = temp_obj._knots[i];
            }
            for (int i = 0; i < temp_obj._knots.size() - knot_span; ++i) {
                surf2_kv[i] = temp_obj._knots[i + knot_span];
            }
            // Add param to the end of surf1_kv and beginning of surf2_kv
            surf1_kv.push_back(param);
            //surf2_kv.insert(surf2_kv.begin(),param);
            for (int j = 0; j <= _degree; ++j) {
                surf2_kv.insert(surf2_kv.begin(), param);
            }

            // Control points for the new curves
            std::vector<vec4> curve1_ctrlpts(temp_obj._control_points.begin(),
                                             temp_obj._control_points.begin() + ks + r);
            std::vector<vec4> curve2_ctrlpts(temp_obj._control_points.begin() + ks + r - 1,
                                             temp_obj._control_points.end());

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
        void evaluate(double t, vec3 &result) const {
            int n = static_cast<int>(_control_points.size()) - 1;
            vec4 res = {0.0, 0.0, 0.0, 0.0};

            // Assume curve_point is already implemented
            curve_point(n, _degree, _knots, _control_points, t, res, _periodic);

            result[0] = res[0];
            result[1] = res[1];
            result[2] = res[2];
        }

        const std::vector<vec4> get_control_points() {
            return _control_points;
        }

        void set_control_points(std::vector<vec4> &cpts) {
            bool change_size = (cpts.size() != _control_points.size());
            _control_points = std::move(cpts);
            if (change_size) {
                if (_periodic) {
                    generate_knots_periodic();
                } else {
                    generate_knots();
                }
            }
        }

        int get_degree() {
            return _degree;
        }

        void set_degree(int val) {
            _degree = val;
            _update_interval();
        }

        const std::vector<double> get_knots() {
            return _knots;
        }

        void set_knots(std::vector<double> &knots) {
            _knots = std::move(knots);
            _update_interval();
        }

        const std::array<double, 2> interval() {
            return _interval;
        }

        AABB &aabb() {
            //rebuildAABB();

            return _aabb;
        }


        // Member variables


        std::vector<vec4> _control_points{}; // Each control point is (x, y, z, weight)
        int _degree = 0;
        std::vector<double> _knots{};
        bool _periodic = false;
        std::vector<vec3> _hull{};
        std::array<double, 2> _interval{0., 1.}; // [min_knot, max_knot]
        AABB _aabb{};
        // Helper methods
        void rebuildAABB() {
            _aabb.min.set(_control_points[0].to_vec3());
            _aabb.max.set(_control_points[0].to_vec3());
            for (auto p: _control_points) {
                _aabb.expand(p.to_vec3());
            }
        }
    };


    class NURBSSurface {
    };
}


#endif //NURBS_H
