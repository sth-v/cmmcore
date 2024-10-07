//
// Created by Andrew Astakhov on 28.09.24.
//

#ifndef NURBS_UTILS_H
#define NURBS_UTILS_H
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <set>
#include <stdexcept>
#include "cmmcore/vec.h"
#include "cmmcore/binom.h"
#include "cmmcore/utils.h"

namespace cmmcore {
    /**
     * @brief Determine the knot span index for a given parameter value `u`.
     *
     * This function finds the knot span index `i` such that the parameter `u`
     * lies within the interval [U[i], U[i+1]] in the knot vector `U`.
     * The knot vector `U` is assumed to be non-decreasing and the parameter
     * `u` is within the range `[U[p], U[n+1]]`.
     *
     * @param n The maximum index of the knot span, typically the number of basis functions minus one.
     * @param p The degree of the B-spline or NURBS.
     * @param u The parameter value for which the span index is to be found.
     * @param U The knot vector, a non-decreasing sequence of real numbers.
     * @param is_periodic Indicates if the curve is periodic.
     * @return The index `i` such that `U[i] <= u < U[i+1]`, where `i` is the knot span.
     * @throws std::invalid_argument If the parameter `u` is outside the bounds of the knot vector `U`.
     *
     * @note
     * The function employs a binary search algorithm to efficiently locate
     * the knot span. It handles special cases where `u` is exactly equal to
     * the last value in `U` and when `u` is outside the range of `U`.
     *
     * @example
     * std::vector<double> U = {0, 0, 0, 0.5, 1, 1, 1};
     * int span = find_span(4, 2, 0.3, U, false);
     * // span should be 2
     */
    inline int find_span(const int n, const int p, double u, const std::vector<double> &U,
                         const bool is_periodic) noexcept {
        //assert(U.size() >= static_cast<size_t>(n + 2));

        const double U_min = U[p];
        const double U_max = U[n + 1];
        double period;

        if (is_periodic) {
            // Wrap u to be within the valid range for periodic and closed curves
            period = U_max - U_min;
            while (u < U_min) {
                u += period;
            }
            while (u > U_max) {
                u -= period;
            }
        } else {
            // Clamp u to be within the valid range for open curves
            if (u >= U[n + 1]) {
                return n;
            } else if (u < U[0]) {
                return p;
            }
        }

        // Handle special case for the upper boundary
        if (u == U[n + 1]) {
            return n;
        }

        // Binary search for the correct knot span
        int low = p;
        int high = n + 1;
        int mid = (low + high) / 2;

        while (u < U[mid] || u >= U[mid + 1]) {
            if (u < U[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low + high) / 2;
        }

        return mid;
    }

    /**
     * @brief Compute the nonvanishing basis functions.
     *
     * @param i The span index.
     * @param u The parameter value.
     * @param p The degree of the B-spline or NURBS.
     * @param U The knot vector.
     * @param N Output vector to store the basis functions.
     */
    inline void basis_funs(const int i, const double u, const int p, const std::vector<double> &U,
                           std::vector<double> &N) noexcept {
        int pp = p + 1;
        N.resize(pp, 0.0);

        std::vector<double> left(pp, 0.0);
        std::vector<double> right(pp, 0.0);
        N[0] = 1.0;

        for (int j = 1; j < pp; ++j) {
            left[j] = u - U[i + 1 - j];
            right[j] = U[i + j] - u;

            double saved = 0.0;
            for (int r = 0; r < j; ++r) {
                double temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
    }

    /**
     * @brief Compute a point on a B-spline curve.
     *
     * @param n The number of basis functions minus one.
     * @param p The degree of the basis functions.
     * @param U The knot vector.
     * @param P The control points of the B-spline curve.
     * @param u The parameter value.
     * @param result Output array to store the computed curve point.
     * @param is_periodic Indicates if the curve is periodic.
     */
    inline void curve_point(const int n, const int p, const std::vector<double> &U,
                            const std::vector<vec4> &P, double u,
                            vec3 &result, const bool is_periodic) noexcept {
        const int pp = p + 1;
        result.set(0.,0.,0.); // Initialize result with 4 zeros

        const int span = find_span(n, p, u, U, is_periodic);
        std::vector<double> N(pp, 0.0);
        basis_funs(span, u, p, U, N);

        double sum_of_weights = 0.0;

        for (int i = 0; i < pp; ++i) {
            const int idx = span - p + i;
            const double b = N[i] * P[idx].w;
            sum_of_weights += b;
            result.x += b * P[idx].x;
            result.y += b * P[idx].y;
            result.z += b * P[idx].z;
        }

        result.x /= sum_of_weights;
        result.y /= sum_of_weights;
        result.z /= sum_of_weights;

    }


    /**
     * @brief Compute the nonzero basis functions and their derivatives for a B-spline.
     *
     * This function calculates the nonzero basis functions and their derivatives
     * for a given parameter value `u` in a B-spline curve. The derivatives are
     * computed up to the `n`-th derivative.
     *
     * @param i The knot span index such that `U[i] <= u < U[i+1]`.
     * @param u The parameter value at which the basis functions and their derivatives are evaluated.
     * @param p The degree of the B-spline basis functions.
     * @param n The number of derivatives to compute (0 for just the basis functions).
     * @param U The knot vector.
     * @param ders Output 2D vector where `ders[k][j]` contains the `k`-th derivative of the `j`-th nonzero basis function at `u`.
     */
    inline void ders_basis_funs(const int i, const double u, const int p, const int n, const std::vector<double> &U,
                                std::vector<vec4> &ders) {
        int pp = p + 1;
        int nn = n + 1;

        // Initialize ders with zeros
        ders.assign(nn, vec4( 0.0,0.0,0.0,.0));

        std::vector<vec4> ndu(pp,  vec4( 0.0,0.0,0.0,.0));
        std::vector<double> left(pp, 0.0);
        std::vector<double> right(pp, 0.0);
        std::vector<vec4> a(2, vec4( 0.0,0.0,0.0,.0));

        ndu[0][0] = 1.0;

        double saved, temp;

        // Compute the basis functions and store in ndu
        for (int j = 1; j < pp; ++j) {
            left[j] = u - U[i + 1 - j];
            right[j] = U[i + j] - u;
            saved = 0.0;

            for (int r = 0; r < j; ++r) {
                ndu[j][r] = right[r + 1] + left[j - r];
                temp = ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        // Load the basis functions
        for (int j = 0; j < pp; ++j) {
            ders[0][j] = ndu[j][p];
        }

        // Compute the derivatives (loop over function index)
        for (int r = 0; r < pp; ++r) {
            int s1 = 0;
            int s2 = 1;
            a[0][0] = 1.0;

            // Loop over derivative order
            for (int k = 1; k <= n; ++k) {
                double d = 0.0;
                int rk = r - k;
                int pk = p - k;

                if (r >= k) {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                    d = a[s2][0] * ndu[rk][pk];
                } else {
                    a[s2][0] = 0.0;
                }

                int j1 = (rk >= -1) ? 1 : -rk;
                int j2 = (r - 1 <= pk) ? k - 1 : p - r;

                for (int j = j1; j <= j2; ++j) {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                    d += a[s2][j] * ndu[rk + j][pk];
                }

                if (r <= pk) {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                    d += a[s2][k] * ndu[r][pk];
                } else {
                    a[s2][k] = 0.0;
                }

                ders[k][r] = d;

                // Swap rows
                std::swap(s1, s2);
            }
        }

        // Multiply through by the correct factors
        int r = p;
        for (int k = 1; k <= n; ++k) {
            for (int j = 0; j < pp; ++j) {
                ders[k][j] *= r;
            }
            r *= (p - k);
        }
    }

    /**
     * @brief Compute the derivatives of a B-spline curve.
     *
     * This function computes the derivatives of a B-spline curve at a given parameter value `u`.
     * The derivatives are computed up to the `d`-th derivative.
     *
     * @param n The number of basis functions minus one.
     * @param p The degree of the basis functions.
     * @param U The knot vector.
     * @param P The control points of the B-spline curve.
     * @param u The parameter value.
     * @param d The number of derivatives to compute.
     * @param CK Output 2D vector to store the computed curve derivatives.
     * @param is_periodic Indicates if the curve is periodic.
     */
    inline void curve_derivs_alg1(const int n, const int p, const std::vector<double> &U,
                                  const std::vector<vec4> &P, double u, int d,
                                  std::vector<vec4> &CK, bool is_periodic) {
        int du = std::min(d, p);

        // Initialize CK with zeros
        CK.assign(du + 1, vec4( 0.0,0.0,0.0,.0));

        int span = find_span(n, p, u, U, is_periodic);

        std::vector<vec4> nders;
        ders_basis_funs(span, u, p, du, U, nders);

        int pp = p + 1;

        // Compute the curve derivatives
        for (int k = 0; k <= du; ++k) {
            for (int j = 0; j < pp; ++j) {
                unsigned long idx = span - p + j;

                // Handle periodicity in control points
                if (is_periodic) {
                    idx = (idx + P.size()) % P.size();
                }

                for (size_t l = 0; l < P[0].size(); ++l) {
                    CK[k][l] += nders[k][j] * P[idx][l];
                }
            }
        }
    }

    /**
     * @brief Compute control points of curve derivatives.
     *
     * This function computes the control points of the derivatives of a B-spline curve
     * up to a specified derivative order.
     *
     * @param p Degree of the basis functions.
     * @param U Knot vector.
     * @param P Control points of the B-spline curve.
     * @param d Number of derivatives to compute.
     * @param r1 Start index for control points.
     * @param r2 End index for control points.
     * @param PK Output 3D vector to store the computed control points of curve derivatives.
     */
    inline void curve_deriv_cpts(const int p, const std::vector<double> &U,
                                 const std::vector<vec4> &P, const int d, const int r1, const int r2,
                                 std::vector<std::vector<vec4> > &PK) {
        int r = r2 - r1;
        int dim = static_cast<int>(P[0].size());
        int pp = p + 1;

        // Initialize PK with zeros
        PK.assign(d + 1, std::vector<vec4>(r + 1, vec4( 0.0,0.0,0.0,.0)));

        // Load the initial control points
        for (int i = 0; i <= r; ++i) {
            PK[0][i] = P[r1 + i];
        }

        // Compute the control points for derivatives
        for (int k = 1; k <= d; ++k) {
            double tmp = p - k + 1;
            for (int i = 0; i <= r - k; ++i) {
                double denom = U[r1 + i + pp] - U[r1 + i + k];
                assert(denom != 0.0 && "Denominator cannot be zero in curve_deriv_cpts.");
                for (int j = 0; j < dim; ++j) {
                    PK[k][i][j] = (tmp * (PK[k - 1][i + 1][j] - PK[k - 1][i][j])) / denom;
                }
            }
        }
    }

    /**
     * @brief Compute the control points of the partial derivatives of a B-spline surface.
     *
     * This function computes the control points of the partial derivatives up to a given order
     * for a surface defined by tensor-product B-splines.
     *
     * @param dim Dimensionality of the control points (e.g., 3 for 3D).
     * @param degree Degrees of the B-spline in the U and V directions.
     *               `degree[0]` is degree in U, `degree[1]` is degree in V.
     * @param kv0 Knot vector for the U direction.
     * @param kv1 Knot vector for the V direction.
     * @param cpts Control points of the B-spline surface.
     *             Indexed as cpts[v][u][dim].
     * @param cpsize Number of control points in U and V directions.
     *               `cpsize[0]` is size in U, `cpsize[1]` is size in V.
     * @param rs Range of indices in the U direction (`rs[0]` to `rs[1]`).
     * @param ss Range of indices in the V direction (`ss[0]` to `ss[1]`).
     * @param deriv_order Order of the highest derivative to compute.
     * @param PKL Output array to store the computed derivatives of the control points.
     *             Indexed as PKL[du+1][dv+1][r+1][s+1][dim].
     */
    inline void surface_deriv_cpts(const int dim, const std::vector<int> &degree,
                                   const std::vector<double> &kv0, const std::vector<double> &kv1,
                                   const std::vector<std::vector<vec4> > &cpts,
                                   const std::vector<int> &cpsize, const std::vector<int> &rs,
                                   const std::vector<int> &ss, const int deriv_order,
                                   std::vector<std::vector<std::vector<std::vector<vec4> > > > &PKL) {
        int du = std::min(degree[0], deriv_order);
        int dv = std::min(degree[1], deriv_order);
        int r = rs[1] - rs[0];
        int s = ss[1] - ss[0];

        // Initialize PKL with zeros
        PKL.assign(du + 1, std::vector<std::vector<std::vector<vec4> > >(
                       dv + 1, std::vector<std::vector<vec4> >(
                           r + 1, std::vector<vec4>(
                               s + 1, vec4( 0.0,0.0,0.0,.0)))));

        std::vector<vec4> temp_cpts(cpsize[0], vec4( 0.0,0.0,0.0,.0));
        std::vector<std::vector<vec4> > PKu;

        // Control points of the U derivatives of every V-curve
        for (int j = ss[0]; j <= ss[1]; ++j) {
            // Extract the j-th row of control points
            for (int i = 0; i < cpsize[0]; ++i) {
                temp_cpts[i] = cpts[j][i];
            }

            // Compute derivative control points in U direction
            curve_deriv_cpts(degree[0], kv0, temp_cpts, du, rs[0], rs[1], PKu);

            // Store the U partial derivatives in PKL
            for (int k = 0; k <= du; ++k) {
                for (int i = 0; i <= r - k; ++i) {
                    PKL[k][0][i][j - ss[0]] = PKu[k][i];
                }
            }
        }

        // Temporary variables
        std::vector<vec4> temp_cpts_v(cpsize[1], vec4( 0.0,0.0,0.0,.0));
        std::vector<std::vector<vec4> > PKuv;

        // Control points of the V derivatives of every U-differentiated V-curve
        for (int k = 0; k <= du; ++k) {
            for (int i = 0; i <= r - k; ++i) {
                int dd = std::min(deriv_order - k, dv);

                // Extract the control points for the current U-curve
                for (int j = 0; j <= s; ++j) {
                    temp_cpts_v[j] = PKL[k][0][i][j];
                }

                // Compute derivative control points in V direction
                curve_deriv_cpts(degree[1], kv1, temp_cpts_v, dd, ss[0], ss[1], PKuv);

                // Store the V partial derivatives in PKL
                for (int l = 1; l <= dd; ++l) {
                    for (int j = 0; j <= s - l; ++j) {
                        PKL[k][l][i][j] = PKuv[l][j];
                    }
                }
            }
        }
    }


    /**
     * @brief Compute a point on a B-spline surface.
     *
     * This function computes a point on a B-spline surface given parameter values (u, v).
     *
     * @param n Number of basis functions in U direction minus one.
     * @param p Degree of the basis functions in U direction.
     * @param U Knot vector in U direction.
     * @param m Number of basis functions in V direction minus one.
     * @param q Degree of the basis functions in V direction.
     * @param V Knot vector in V direction.
     * @param Pw Control points of the B-spline surface (with weights), indexed as Pw[u][v][dim].
     * @param u Parameter value in U direction.
     * @param v Parameter value in V direction.
     * @param periodic_u Indicates if the surface is periodic in U direction.
     * @param periodic_v Indicates if the surface is periodic in V direction.
     * @param result Output vector to store the computed surface point (dimension 3).
     */
    inline void surface_point(const int n, const int p, const std::vector<double> &U, const int m, const int q,
                              const std::vector<double> &V,
                              const std::vector<std::vector<vec4> > &Pw, double u,
                              const double v, const bool periodic_u, const bool periodic_v,
                              vec3 &result) {
        int uspan = find_span(n, p, u, U, periodic_u);
        int vspan = find_span(m, q, v, V, periodic_v);

        std::vector<double> Nu(p + 1, 0.0);
        std::vector<double> Nv(q + 1, 0.0);

        basis_funs(uspan, u, p, U, Nu);
        basis_funs(vspan, v, q, V, Nv);

        const int dim = static_cast<int>(Pw[0][0].size());
        std::vector<vec4> temp(q + 1, vec4( 0.0,0.0,0.0,.0));

        // Compute the temporary points
        for (int l = 0; l <= q; ++l) {
            for (int k = 0; k <= p; ++k) {
                unsigned long u_idx = uspan - p + k;
                unsigned long v_idx = vspan - q + l;

                // Handle periodicity in control points
                if (periodic_u) {
                    u_idx = (u_idx + Pw.size()) % Pw.size();
                }
                if (periodic_v) {
                    v_idx = (v_idx + Pw[0].size()) % Pw[0].size();
                }

                for (int d = 0; d < dim; ++d) {
                    temp[l][d] += Nu[k] * Pw[u_idx][v_idx][d];
                }
            }
        }

        // Compute the surface point
        std::vector<double> Sw(dim, 0.0);
        for (int l = 0; l <= q; ++l) {
            for (int d = 0; d < dim; ++d) {
                Sw[d] += Nv[l] * temp[l][d];
            }
        }

        // Dehomogenize the point
        double w = Sw[dim - 1];
        assert(w != 0.0 && "Weight w cannot be zero in surface_point.");


        for (int d = 0; d < dim - 1; ++d) {
            result[d] = Sw[d] / w;
        }
    }


    /**
     * @brief Computes the derivatives of a rational B-spline surface.
     *
     * This function calculates the derivatives of a rational B-spline surface based on the given
     * weighted surface derivative values (`SKLw`). The results are stored in the `SKL` array.
     *
     * @param SKLw Input array containing the weighted surface derivatives.
     *             Indexed as SKLw[k][l][dim].
     * @param deriv_order Highest order of derivative to compute.
     * @param SKL Output array to store the computed derivatives of the surface.
     *            Indexed as SKL[k][l][dim - 1].
     */
    inline void rat_surface_derivs(const std::vector<std::vector<vec4> > &SKLw,
                                   const int deriv_order, std::vector<std::vector<vec4> > &SKL) {
        const auto dim = SKLw[0][0].size();
        const auto dm = dim - 1;
        const auto do_order = deriv_order + 1;

        // Initialize SKL with zeros
        SKL.assign(do_order, std::vector<vec4>(do_order, vec4( 0.0,0.0,0.0,.0)));

        // Temporary variables
        vec4 v(0.0, 0.0,0.0,0.0);
        vec3 v2(0.0, 0.0,0.0);


        for (int k = 0; k < do_order; ++k) {
            for (int l = 0; l < do_order; ++l) {
                v = SKLw[k][l];

                // Subtract contributions from lower-order derivatives
                for (int j = 1; j <= l; ++j) {
                    const double bin_lj = binomial_coefficient(l, j);
                    for (int d = 0; d < dm; ++d) {
                        v[d] -= bin_lj * SKLw[0][j][dm] * SKL[k][l - j][d];
                    }
                }

                for (int i = 1; i <= k; ++i) {
                    const double bin_ki = binomial_coefficient(k, i);

                    for (int d = 0; d < dm; ++d) {
                        v[d] -= bin_ki * SKLw[i][0][dm] * SKL[k - i][l][d];
                    }

                    // Reset v2
                    v2.set(0,0,0);


                    for (int j = 1; j <= l; ++j) {
                        double bin_lj = binomial_coefficient(l, j);
                        for (int d = 0; d < dm; ++d) {
                            v2[d] += bin_lj * SKLw[i][j][dm] * SKL[k - i][l - j][d];
                        }
                    }

                    for (int d = 0; d < dm; ++d) {
                        v[d] -= bin_ki * v2[d];
                    }
                }

                // Divide by weight
                double w = SKLw[0][0][dm];
                assert(w != 0.0 && "Weight w cannot be zero in rat_surface_derivs.");

                for (int d = 0; d < dm; ++d) {
                    SKL[k][l][d] = v[d] / w;
                }
            }
        }
    }


    /**
     * @brief Finds the multiplicity of a knot in a knot vector.
     *
     * @param knot The knot value to find.
     * @param knot_vector The knot vector.
     * @return The multiplicity of the knot.
     */
    inline int find_multiplicity(const double knot, const std::vector<double> &knot_vector) noexcept {
        int mult = 0;
        for (const auto &kv: knot_vector) {
            double difference = knot - kv;
            if (std::fabs(difference) <= 1e-12) {
                mult += 1;
            }
        }
        return mult;
    }

    /**
 * @brief Computes the alpha value for knot insertion.
 *
 * @param u The knot value to insert.
 * @param knots The knot vector.
 * @param span The knot span index.
 * @param idx The current index.
 * @param leg The leg value.
 * @return The alpha value.
 */
    inline double knot_insertion_alpha(const double u, const std::vector<double> &knots, const int span, const int idx,
                                       const int leg) noexcept {
        return (u - knots[leg + idx]) / (knots[idx + span + 1] - knots[leg + idx]);
    }

    /**
     * @brief Computes the alpha_i value for knot removal.
     *
     * @param u The knot value.
     * @param degree The degree of the B-spline.
     * @param knots The knot vector.
     * @param num The current iteration number.
     * @param idx The index i.
     * @return The alpha_i value.
     */
    inline double knot_removal_alpha_i(const double u, const int degree, const std::vector<double> &knots,
                                       const int num, const int idx) noexcept {
        return (u - knots[idx]) / (knots[idx + degree + 1 + num] - knots[idx]);
    }


    /**
     * @brief Computes the alpha_j value for knot removal.
     *
     * @param u The knot value.
     * @param degree The degree of the B-spline.
     * @param knots The knot vector.
     * @param num The current iteration number.
     * @param idx The index j.
     * @return The alpha_j value.
     */
    inline double knot_removal_alpha_j(const double u, const int degree, const std::vector<double> &knots,
                                       const int num, const int idx) noexcept {
        return (u - knots[idx - num]) / (knots[idx + degree + 1] - knots[idx - num]);
    }


    /**
     * @brief Performs knot insertion in a B-spline curve.
     *
     * @param degree Degree of the B-spline.
     * @param knots Knot vector.
     * @param ctrlpts Control points of the B-spline curve.
     * @param u Knot value to insert.
     * @param num Number of times to insert the knot.
     * @param s Current multiplicity of the knot.
     * @param span Knot span index.
     * @param is_periodic Indicates if the curve is periodic.
     * @return New control points after knot insertion.
     */
    inline std::vector<vec4> knot_insertion(const int degree,
                                                            const std::vector<double> &knots,
                                                            const std::vector<vec4> &ctrlpts,
                                                            const double u, const int num, const int s, const int span
    ) {

        int n = static_cast<int>(ctrlpts.size());//static_cast<int>(ctrlpts.size()) - 1;
        int nq = n + num; //int nq = n + num + 1;

        int dim = static_cast<int>(ctrlpts[0].size());

        std::vector<vec4> result(nq, vec4(0.,0.,0.,0.));
        std::vector<vec4> temp(degree + 1, vec4(0.,0.,0.,0.));

        // Copy unaffected control points
        for (int i = 0; i <= span - degree; ++i) {
            result[i] = ctrlpts[i];
        }
        for (int i = span - s; i < n; ++i) {//for (int i = span - s + 1; i <= n; ++i) {
            result[i + num] = ctrlpts[i];
        }

        // Copy affected control points
        for (int i = 0; i <= degree - s; ++i) {
            temp[i] = ctrlpts[span - degree + i];
        }

        // Knot insertion algorithm
        for (int j = 1; j <= num; ++j) {
            int L = span - degree + j;
            for (int i = 0; i <= degree - j - s; ++i) {
                double alpha = knot_insertion_alpha(u, knots, span, i, L);
                for (int idx = 0; idx < dim; ++idx) {
                    temp[i][idx] = alpha * temp[i + 1][idx] + (1.0 - alpha) * temp[i][idx];
                }
            }
            //memcpy(&result[L][ 0], &temp[0][0], sizeof(double) * dim);
            //memcpy(&result[span + num - j - s][ 0], &temp[(degree - j - s)][0], sizeof(double) * dim);
            result[L] = temp[0];

            result[span + num - j - s] = temp[degree - j - s];


        }

        // Copy remaining affected control points
        int L = span - degree + num;
        for (int i = L + 1; i < span - s + 1; ++i) {
            result[i] = temp[i - L];
        }

        return result;
    }
    inline void knot_insertion(const int degree,
                                                            const std::vector<double> &knots,
                                                            const std::vector<vec4> &ctrlpts,
                                                            const double u, const int num, const int s, const int span,std::vector<vec4>& result
    ) {

        int n = static_cast<int>(ctrlpts.size());//static_cast<int>(ctrlpts.size()) - 1;
        //int nq = n + num; //int nq = n + num + 1;

        int dim = static_cast<int>(ctrlpts[0].size());


        std::vector<vec4> temp(degree + 1, vec4(0.,0.,0.,0.));

        // Copy unaffected control points
        for (int i = 0; i <= span - degree; ++i) {
            result[i] = ctrlpts[i];
        }
        for (int i = span - s; i < n; ++i) {//for (int i = span - s + 1; i <= n; ++i) {
            result[i + num] = ctrlpts[i];
        }

        // Copy affected control points
        for (int i = 0; i <= degree - s; ++i) {
            temp[i] = ctrlpts[span - degree + i];
        }

        // Knot insertion algorithm
        for (int j = 1; j <= num; ++j) {
            int L = span - degree + j;
            for (int i = 0; i <= degree - j - s; ++i) {
                double alpha = knot_insertion_alpha(u, knots, span, i, L);
                for (int idx = 0; idx < dim; ++idx) {
                    temp[i][idx] = alpha * temp[i + 1][idx] + (1.0 - alpha) * temp[i][idx];
                }
            }
            //memcpy(&result[L][ 0], &temp[0][0], sizeof(double) * dim);
            //memcpy(&result[span + num - j - s][ 0], &temp[(degree - j - s)][0], sizeof(double) * dim);
            result[L] = temp[0];

            result[span + num - j - s] = temp[degree - j - s];


        }

        // Copy remaining affected control points
        int L = span - degree + num;
        for (int i = L + 1; i < span - s + 1; ++i) {
            result[i] = temp[i - L];
        }


    }

    /**
 * @brief Inserts knots into a knot vector.
 *
 * @param knots Original knot vector.
 * @param u Knot value to insert.
 * @param span Knot span index.
 * @param r Number of times to insert the knot.
 * @return New knot vector after insertion.
 */
    inline std::vector<double> knot_insertion_kv(const std::vector<double> &knots, const double u, int span, int r) {
        //auto kv_size =knots.size();
        std::vector<double> kv_updated=knots;



        for (int i = 0; i < r; ++i) {
            kv_updated.insert(kv_updated.begin()+span+1, u);
        }

        return kv_updated;
    }


    /**
     * @brief Calculates the Euclidean distance between two points in N-dimensional space.
     *
     * @param a Pointer to the first point.
     * @param b Pointer to the second point.
     * @param dim Dimension of the points.
     * @return The Euclidean distance between the points.
     */
    inline double point_distance(const double *a, const double *b, const int dim) noexcept {
        double temp = 0.0;
        for (int i = 0; i < dim; ++i) {
            temp += std::pow(a[i] - b[i], 2);
        }
        return std::sqrt(temp);
    }


    /**
     * @brief Performs knot removal from a B-spline curve.
     *
     * @param degree Degree of the B-spline.
     * @param knots Knot vector.
     * @param ctrlpts Control points of the B-spline curve.
     * @param u Knot value to remove.
     * @param tol Tolerance for point distance comparison.
     * @param num Number of times to remove the knot.
     * @param is_periodic Indicates if the curve is periodic.
     * @return New control points after knot removal.
     */
    inline std::vector<vec4> knot_removal(int degree, const std::vector<double> &knots,
                                                          const std::vector<vec4> &ctrlpts,
                                                          double u, double tol = 1e-4, int num = 1,
                                                          bool is_periodic = false) {
        int s = find_multiplicity(u, knots);
        int n = static_cast<int>(ctrlpts.size()) - 1;
        int r = find_span(n, degree, u, knots, is_periodic);

        int first = r - degree;
        int last = r - s;
        int dim = static_cast<int>(ctrlpts[0].size());

        std::vector<vec4> ctrlpts_new = ctrlpts;
        std::vector<vec4> temp(2 * degree + 1,vec4(0.,0.,0.,0.));

        for (int t = 0; t < num; ++t) {
            temp[0] = ctrlpts[first - 1];
            temp[last - first + 2] = ctrlpts[last + 1];
            int i = first;
            int j = last;
            int ii = 1;
            int jj = last - first + 1;
            bool remflag = false;

            while (j - i >= t) {
                double alpha_i = knot_removal_alpha_i(u, degree, knots, t, i);
                double alpha_j = knot_removal_alpha_j(u, degree, knots, t, j);
                for (int k = 0; k < dim; ++k) {
                    temp[ii][k] = (ctrlpts[i][k] - (1.0 - alpha_i) * temp[ii - 1][k]) / alpha_i;
                    temp[jj][k] = (ctrlpts[j][k] - alpha_j * temp[jj + 1][k]) / (1.0 - alpha_j);
                }
                i += 1;
                j -= 1;
                ii += 1;
                jj -= 1;
            }

            vec4 ptn(0.0,0.,0.,0.);
            if (j - i < t) {
                if (point_distance(&temp[ii - 1][0], &temp[jj + 1][0], dim) <= tol) {
                    remflag = true;
                }
            } else {
                double alpha_i = knot_removal_alpha_i(u, degree, knots, t, i);
                for (int k = 0; k < dim; ++k) {
                    ptn[k] = alpha_i * temp[ii + t + 1][k] + (1.0 - alpha_i) * temp[ii - 1][k];
                }
                if (point_distance(&ctrlpts[i].x, &ptn.x, dim) <= tol) {
                    remflag = true;
                }
            }

            if (remflag) {
                i = first;
                j = last;
                while (j - i > t) {
                    ctrlpts_new[i] = temp[i - first + 1];
                    ctrlpts_new[j] = temp[j - first + 1];
                    i += 1;
                    j -= 1;
                }
            }

            first -= 1;
            last += 1;
        }

        int t_total = num;
        int new_size = static_cast<int>(ctrlpts_new.size()) - t_total;

        ctrlpts_new.resize(new_size);

        return ctrlpts_new;
    }


    /**
     * @brief Performs knot refinement on a B-spline curve.
     *
     * @param degree Degree of the B-spline.
     * @param knots Knot vector.
     * @param ctrlpts Control points of the B-spline curve.
     * @param knot_list List of knots in the parameter range.
     * @param add_knot_list Additional knots to be added.
     * @param density Refinement density.
     * @param is_periodic Indicates if the curve is periodic.
     * @return A pair containing new control points and new knot vector after refinement.
     */
    inline std::pair<std::vector<vec4>, std::vector<double> > knot_refinement(
        int degree, const std::vector<double> &knots,
        const std::vector<vec4> &ctrlpts,
        const std::vector<double> &knot_list = {},
        const std::vector<double> &add_knot_list = {},
        int density = 1, bool is_periodic = false) {
        int n = static_cast<int>(ctrlpts.size()) - 1;
        int m = n + degree + 1;
        int dim = static_cast<int>(ctrlpts[0].size());

        // Combine knot_list and add_knot_list
        std::vector<double> new_knots = knot_list.empty()
                                            ? std::vector<double>(knots.begin() + degree, knots.end() - degree)
                                            : knot_list;

        if (!add_knot_list.empty()) {
            new_knots.insert(new_knots.end(), add_knot_list.begin(), add_knot_list.end());
        }

        // Remove duplicates and sort the knot list
        std::set<double> knot_set(new_knots.begin(), new_knots.end());
        new_knots.assign(knot_set.begin(), knot_set.end());

        // Apply density refinement
        for (int d = 0; d < density - 1; ++d) {
            std::vector<double> rknots;
            for (size_t i = 0; i < new_knots.size() - 1; ++i) {
                rknots.push_back(new_knots[i]);
                rknots.push_back(new_knots[i] + (new_knots[i + 1] - new_knots[i]) / 2.0);
            }
            rknots.push_back(new_knots.back());
            new_knots = rknots;
        }

        // Determine the knots to be inserted
        std::vector<double> X;
        for (const auto &mk: new_knots) {
            int s = find_multiplicity(mk, knots);
            int r = degree - s;
            for (int i = 0; i < r; ++i) {
                X.push_back(mk);
            }
        }

        if (X.empty()) {
            throw std::runtime_error("Cannot refine knot vector on this parametric dimension");
        }

        int r_val = static_cast<int>(X.size()) - 1;
        int a = find_span(n, degree, X.front(), knots, is_periodic);
        int b = find_span(n, degree, X.back(), knots, is_periodic) + 1;

        std::vector<vec4> new_ctrlpts(n + r_val + 2, vec4(0.,0.,0.,0.));
        std::vector<double> new_kv(m + r_val + 2, 0.0);

        // Copy unaffected control points and knot vector
        for (int j = 0; j <= a - degree; ++j) {
            new_ctrlpts[j] = ctrlpts[j];
        }
        for (int j = b - 1; j <= n; ++j) {
            new_ctrlpts[j + r_val + 1] = ctrlpts[j];
        }
        for (int j = 0; j <= a; ++j) {
            new_kv[j] = knots[j];
        }
        for (int j = b + degree; j <= m; ++j) {
            new_kv[j + r_val + 1] = knots[j];
        }

        int i = b + degree - 1;
        int k = b + degree + r_val;
        int j = r_val;

        while (j >= 0) {
            while (X[j] <= knots[i] && i > a) {
                new_ctrlpts[k - degree - 1] = ctrlpts[i - degree - 1];
                new_kv[k] = knots[i];
                k -= 1;
                i -= 1;
            }
            new_ctrlpts[k - degree - 1] = new_ctrlpts[k - degree];
            for (int l = 1; l <= degree; ++l) {
                int idx = k - degree + l;
                double alpha = new_kv[k + l] - X[j];
                if (std::fabs(alpha) < 1e-8) {
                    new_ctrlpts[idx - 1] = new_ctrlpts[idx];
                } else {
                    alpha = alpha / (new_kv[k + l] - knots[i - degree + l]);
                    for (int idx2 = 0; idx2 < dim; ++idx2) {
                        new_ctrlpts[idx - 1][idx2] =
                                alpha * new_ctrlpts[idx - 1][idx2] + (1.0 - alpha) * new_ctrlpts[idx][idx2];
                    }
                }
            }
            new_kv[k] = X[j];
            k -= 1;
            j -= 1;
        }

        return {new_ctrlpts, new_kv};
    }
}

#endif //NURBS_UTILS_H
