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
#include <cassert>
#include <unordered_set>

#ifdef CYTHON_ABI
#include "nurbs_utils.h"
#include "utils.h"
#include "aabb.h"
#else
#include "cmmcore/nurbs_utils.h"
#include "cmmcore/aabb.h"
#include "cmmcore/utils.h"
#endif


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
        bool is_periodic() const {
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

            curve_derivs_alg1(n, _degree, knots, control_points, t, 1,     CK, _periodic);
                auto der=CK[1];

            //printf("(%f,%f,%f,%f)\n",der.x,der.y,der.z,der.w);
            result.x= der.x;
            result.y= der.y;
            result.z= der.z;

        }


        void derivative(const double t, std::vector<vec3> &result, const size_t d) const {
            int n = static_cast<int>(control_points.size()) - 1;
            result.resize(d);

            std::vector<vec4> CK(result.size());
            curve_derivs_alg1(n, _degree, knots, control_points, t, d,     CK, _periodic);
            for (int i = 0; i < d; ++i) {
                auto &der=CK[i];
                auto& target=result[i];
                target.x= der.x;
                target.y= der.y;
                target.z= der.z;


            }

        }
        void derivative(const double t, std::vector<vec3> &result) const {
                derivative(t,result,result.size());
        }
        void derivatives(const double t, vec3 &d0, vec3 &d1, vec3 &d2 ) const {
            vec3 temp;
            evaluate(t,d0);
            evaluate(t-1e-5,temp);
            evaluate(t+1e-5,d2);
            d1.set(d2-temp);
            d1/=2e-5;
            d2-=(2*d0);
            d2+=temp;
            d2/=1e-10;


        }
        void derivative_fdm(const double t, vec3 &result) const {
            vec3 temp;
            evaluate(t-1e-5,temp);
            evaluate(t+1e-5,result);
            result-=temp;
            result/=2e-5;
        }
        /**
         ..math:
         \NURBSCurve
         \frac{d^2 f(t)}{dt^2} \approx \frac{f(t + h) - 2f(t) + f(t - h)}{h^2}
         */
        void second_derivative_fdm(const double t, vec3 &result) const {
            vec3 temp;
            evaluate(t-1e-5,temp);
            evaluate(t+1e-5,result);
            result-=temp;
            result/=2e-5;


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

        std::vector<vec4> get_control_points() const {
            return control_points;
        }

        void set_control_points(std::vector<vec4> &cpts) {
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

        int get_degree() {
            return _degree;
        }

        void set_degree(int val) {
            _degree = val;
            update_interval();
        }

        const std::vector<double> get_knots() {
            return knots;
        }

        void set_knots(std::vector<double> &knts) {
            knots = std::move(knts);
            update_interval();
        }

        const std::array<double, 2> interval() {
            return _interval;
        }

        AABB &aabb() {
            //rebuildAABB();

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
    };

    // External function declarations

    inline void flipControlPointsU(const std::vector<std::vector<vec4>> &cpts,std::vector<std::vector<vec4>> &result )

    {
        size_t size_u=cpts.size();
        size_t size_v=cpts[0].size();
        result.clear();


        for (size_t i = 0; i <  size_v; i++) {
            std::vector<vec4> new_row{};
            for (size_t j = 0; j < size_u; j++) {
                new_row.push_back(cpts[j][i]);


            }
            result.push_back(new_row);

        }

    }
    class NURBSSurface {
    public:
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
        void insert_knot(const double t, const int direction, const int count ) {


            if (direction == 0) {
                int s = find_multiplicity(t,  _knots_u);
                int span = find_span(_control_points.size() - 1, _degree[0], t, _knots_u, false);

                auto kv = knot_insertion_kv(_knots_u, t, span, count);
                std::vector<std::vector<vec4> >  cpts_tmp;
                std::vector<vec4>cpts=control_points_flat();

                for (size_t v = 0; v < _size[1]; ++v)
                {
                    std::vector<vec4> ccu;
                    for (size_t u = 0; u < _size[0]; ++u)
                    {
                       ccu.push_back(cpts[u*_size[1] +v]);

                    }


                    cpts_tmp.push_back(knot_insertion(_degree[0], _knots_u, ccu, t,count, s, span));

                }




                flipControlPointsU(cpts_tmp, _control_points);

                _size[0]+=count;

                _knots_u=kv;



            } else {
                insert_knot_v(t, count);
            }}
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

        std::pair<NURBSSurface, NURBSSurface> split_surface_u(double param, double tol = 1e-7) {
            if (param <= _interval[0][0] || param >= _interval[0][1] || std::fabs(param - _interval[0][0]) <= tol ||
                std::fabs(param - _interval[0][1]) <= tol) {
                throw std::invalid_argument("Cannot split from the domain edge");
            }
            int ks = find_span(static_cast<int>(_size[0]-1), _degree[0], param, _knots_u, false) - _degree[0] + 1;
            int s = find_multiplicity(param, _knots_u);

            int r = _degree[0] - s;
            NURBSSurface temp_obj = *this;

            temp_obj.insert_knot_u(param, r);
            std::vector<double> temp_knots_u =temp_obj._knots_u;
            std::vector<std::vector<vec4>> tcpts=temp_obj._control_points;

            int knot_span = find_span(static_cast<int>(temp_obj._size[0]-1), _degree[0], param,temp_knots_u,
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

        std::pair<NURBSSurface, NURBSSurface> split_surface_v(double param, double tol = 1e-7) {
            if (param <= _interval[1][0] || param >= _interval[1][1] || std::fabs(param - _interval[1][0]) <= tol ||
                std::fabs(param - _interval[1][1]) <= tol) {
                throw std::invalid_argument("Cannot split from the domain edge");
            }
            int ks = find_span(static_cast<int>(_size[1]-1), _degree[1], param, _knots_v, false) - _degree[1] + 1;
            int s = find_multiplicity(param, _knots_v);
            int r = _degree[1] - s;
            NURBSSurface temp_obj = *this;
            temp_obj.insert_knot_v(param, r);
            int knot_span = find_span(static_cast<int>(temp_obj._size[1]-1), _degree[1], param, temp_obj._knots_v,
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
        void control_points_flat3d(std::vector<vec3>& control_points_flt ) const {
            control_points_flt.resize(_size[0] * _size[1]);
            for (std::size_t i = 0; i < _size[0]; ++i) {
                for (std::size_t j = 0; j < _size[1]; ++j) {
                    const double w = _control_points[i][j].w;
                    control_points_flt[i * _size[1] + j].set(_control_points[i][j].x / w, _control_points[i][j].y / w,
                                                             _control_points[i][j].z / w);
                }
            };

        }
        std::vector<vec3> control_points_flat3d() const
        {
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


        void subdivide(
            NURBSSurface &s11,
            NURBSSurface &s12,
            NURBSSurface &s21,
            NURBSSurface &s22,
            double umid,
            double vmid
            ) {
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

            ) {
            auto umid = 0.5 * (_interval[0][1] + _interval[0][0]);
            auto vmid = 0.5 * (_interval[1][1] + _interval[1][0]);
            subdivide(s11,s12,s21,s22,umid,vmid);
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

    inline bool contains(std::unordered_set<double> &s, double val) {
        return std::find(s.begin(), s.end(), val) != s.end();
    }

    inline void uniqueKnots(std::vector<double> &knots, std::vector<double> &uniqueKnots) {
        std::unordered_set<double> ks;
        for (std::size_t k = 0; k < knots.size(); ++k) {
            if (!contains(ks, knots[k])) {
                ks.insert(knots[k]);
                uniqueKnots.push_back(knots[k]);
            }
        }
    }

    inline void decomposeDirection(NURBSSurface &surf, std::vector<NURBSSurface> &bezierSurfaces, int direction) {
        //std::vector<double>::iterator knots_it =iterate_unique_knots(_knots_u,_degree[0]);
        assert(direction==0||direction==1);
        NURBSSurface temp = surf;
        std::vector<double> uknots;
        temp.update_interval();
        if (direction == 0) {
            uniqueKnots(temp._knots_u, uknots);

            for (auto &knot: uknots) {
                if (knot != uknots[0] && knot != (uknots[uknots.size() - 1])) {
                    auto [first, second] = temp.split_surface_u(knot);

                    bezierSurfaces.push_back(first);
                    temp = second;
                }
            }
            bezierSurfaces.push_back(temp);
        } else if (direction == 1) {
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

    inline void decompose(NURBSSurface &surf, std::vector<NURBSSurface> &bezierSurfaces) {
        std::vector<NURBSSurface> bezierSurfaces1;
        decomposeDirection(surf, bezierSurfaces1, 0);
        for (auto &item: bezierSurfaces1) {
            decomposeDirection(item, bezierSurfaces, 1);
        }
    }
}


#endif //NURBS_H
