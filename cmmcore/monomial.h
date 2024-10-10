//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef MONOMIAL_H
#define MONOMIAL_H

#ifdef CYTHON_ABI
#include "vec.h"
#include "matrix.h"
#include "binom.h"
#else
#include "cmmcore/vec.h"
#include "cmmcore/matrix.h"
#include "cmmcore/binom.h"
#include "cmmcore/nurbs.h"
#endif

namespace cmmcore {
    enum SurfaceParameter {
        U,
        V
    };


    inline Matrix bpmat(int n) {
        Matrix result(n + 1, n + 1);
        for (int j = 0; j <= n; ++j) {
            for (int k = 0; k <= j; ++k) {
                result(j, k) = std::pow(-1, j - k) * binomial_coefficient(n, j) * binomial_coefficient(j, k);
            }
        }
        return result;
    }

    template<typename T>
    using Tensor2D = std::vector<std::vector<T> >;
#define CMMCORE_CREATE_TENSOR2D(n,m,typ,...) std::vector(n, std::vector<typ>(m, {__VA_ARGS__}))
#define CMMCORE_RESIZE_TENSOR2D(var,n,m,typ,...) var.resize(n, std::vector<typ>(m, {__VA_ARGS__}))

    template<typename T>
    inline void getVectorTensor2DComponent(Tensor2D<T> &arr, int component, Matrix &result) {
        result.rows = arr.size();
        result.cols = arr[0].size();
        result.data.resize(result.cols * result.rows);
        for (int i = 0; i < arr.size(); ++i) {
            for (int j = 0; j < arr[i].size(); ++j) {
                result(i, j) = arr[i][j][component];
            }
        }
    }

    template<typename T>
    inline void setVectorTensor2DComponent(Tensor2D<T> &arr, int component, const Matrix &value) {
        for (int i = 0; i < arr.size(); ++i) {
            for (int j = 0; j < arr[i].size(); ++j) {
                arr[i][j][component] = value(i, j);
            }
        }
    }

    inline void compute_partial_derivative(
        const Tensor2D<vec3> &coeffs,
        const SurfaceParameter variable,
        Tensor2D<vec3> &result) {
        int n = coeffs.size();
        int m = coeffs[0].size();


        if (variable == U) {
            CMMCORE_RESIZE_TENSOR2D(result, n-1, m, vec3, 0, 0, 0);
            for (int i = 1; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    result[i - 1][j][0] = i * coeffs[i][j][0];
                    result[i - 1][j][1] = i * coeffs[i][j][1];
                    result[i - 1][j][2] = i * coeffs[i][j][2];
                }
            }
        } else if (variable == V) {
            CMMCORE_RESIZE_TENSOR2D(result, n, m-1, vec3, 0, 0, 0);

            for (int i = 0; i < n; ++i) {
                for (int j = 1; j < m; ++j) {
                    result[i][j - 1][0] = j * coeffs[i][j][0];
                    result[i][j - 1][1] = j * coeffs[i][j][1];
                    result[i][j - 1][2] = j * coeffs[i][j][2];
                }
            }
        }
    }

    inline void cross(const Tensor2D<vec3> &a_coeffs, const Tensor2D<vec3> &b_coeffs, Tensor2D<vec3> &result) {
        int n1 = a_coeffs.size();
        int m1 = a_coeffs[0].size();
        int n2 = b_coeffs.size();
        int m2 = b_coeffs[0].size();
        int n = n1 + n2 - 1;
        int m = m1 + m2 - 1;

        // Initialize the result array with zeros
        CMMCORE_RESIZE_TENSOR2D(result, n, m, vec3, 0, 0, 0);


        // Compute the cross product
        for (int i1 = 0; i1 < n1; ++i1) {
            for (int j1 = 0; j1 < m1; ++j1) {
                for (int i2 = 0; i2 < n2; ++i2) {
                    for (int j2 = 0; j2 < m2; ++j2) {
                        int i = i1 + i2;
                        int j = j1 + j2;
                        result[i][j][0] += (a_coeffs[i1][j1][1] * b_coeffs[i2][j2][2] - a_coeffs[i1][j1][2] * b_coeffs[
                                                i2][j2][1]);
                        result[i][j][1] += (a_coeffs[i1][j1][2] * b_coeffs[i2][j2][0] - a_coeffs[i1][j1][0] * b_coeffs[
                                                i2][j2][2]);
                        result[i][j][2] += (a_coeffs[i1][j1][0] * b_coeffs[i2][j2][1] - a_coeffs[i1][j1][1] * b_coeffs[
                                                i2][j2][0]);
                    }
                }
            }
        }
    }

    class Monomial2D {
    public:
        int n=0;
        int m=0;
        Tensor2D<vec3> coefficients{};
        Matrix Mu, Mv, Mu_inv, Mv_inv;

        Monomial2D() = default;

        Monomial2D(int _n, int _m, const Tensor2D<vec3> &coeffs): n(_n),m(_m),coefficients(std::move(coeffs)),
                                                                  Mu(bpmat(n - 1)),
                                                                  Mv(bpmat(m - 1)),
                                                                  Mu_inv(std::move(Mu.inverse())),
                                                                  Mv_inv(std::move(Mv.inverse())) {
        }


        Monomial2D(int _n, int _m): n(_n), m(_m), coefficients(std::move(CMMCORE_CREATE_TENSOR2D(_n, _m, vec3, 0, 0, 0))),
                                    Mu(bpmat(n - 1)), Mv(bpmat(m - 1)), Mu_inv(std::move(Mu.inverse())),
                                    Mv_inv(std::move(Mv.inverse())) {
        }

         Monomial2D(const Tensor2D<vec3> &coeffs): n(coeffs.size()), m(coeffs[0].size()), coefficients(std::move(coeffs)),
                                            Mu(bpmat(n - 1)), Mv(bpmat(m - 1)), Mu_inv(std::move(Mu.inverse())),
                                            Mv_inv(std::move(Mv.inverse())) {
        }

        Monomial2D(const Monomial2D &m): n(m.n), m(m.m), coefficients(m.coefficients), Mu(m.Mu), Mv(m.Mv),
                                         Mu_inv(m.Mu_inv), Mv_inv(m.Mv_inv) {
        }

        Monomial2D(const NURBSSurface &surf): n(surf._size[0]), m(surf._size[1]),
                                              coefficients(std::move(CMMCORE_CREATE_TENSOR2D(n, m, vec3, 0, 0, 0))),
                                              Mu(bpmat(n - 1)), Mv(bpmat(m - 1)), Mu_inv(std::move(Mu.inverse())),
                                              Mv_inv(std::move(Mv.inverse())) {
            Tensor2D<vec3> cpts = surf.control_points3d();
            Matrix MvT = Mv.transpose();
            for (int i = 0; i < 3; ++i) {
                Matrix cfs(n, m);
                getVectorTensor2DComponent<vec3>(cpts, i, cfs);
                setVectorTensor2DComponent<vec3>(coefficients, i, Mu * cfs * MvT);
                // Mu @ control_points[:, :, d] @ Mv.T
            }
        }

        void set(Tensor2D<vec3> &coeffs) {
            n = coeffs.size();
            m = coeffs[0].size();
            coefficients = std::move(coeffs);
            Mu = bpmat(n - 1);
            Mv = bpmat(m - 1);
            Mu_inv = std::move(Mu.inverse());
            Mv_inv = std::move(Mv.inverse());
        }

        void to_bezier(NURBSSurface &surf) {
            /*
            *    n, m, dim = monomial_coeffs.shape
            #print(n - 1)
            #print(m-1)
            Mu_inv = np.linalg.inv(bpmat(n - 1, method=bmethod))
            Mv_inv = np.linalg.inv(bpmat(m - 1, method=bmethod))

            control_points = np.zeros((n, m, dim))
            for d in range(dim):
            control_points[:, :, d] = Mu_inv @ monomial_coeffs[:, :, d] @ Mv_inv.T*/
            CMMCORE_RESIZE_TENSOR2D(surf._control_points, n, m, vec4, 0, 0, 0, 1);
            Matrix Mv_inv_T = Mv_inv.transpose();
            for (int i = 0; i < 3; ++i) {
                Matrix cfs(n, m);
                getVectorTensor2DComponent<vec3>(coefficients, i, cfs);
                setVectorTensor2DComponent<vec4>(surf._control_points, i, Mu_inv * cfs * Mv_inv_T);
            }
            surf._degree[0] = n - 1;
            surf._degree[1] = m - 1;
            surf._size[0] = n;
            surf._size[1] = m;
            surf.generate_knots_u();
            surf.generate_knots_v();
            surf._update_interval();
        }

        void computePartialDerivative(const SurfaceParameter variable, Monomial2D &result) const {
            Tensor2D<vec3> coeffs;
            compute_partial_derivative(this->coefficients, variable, coeffs);
            result.set(coeffs);
        }

        void computeNormal(Monomial2D &result) const {
            Tensor2D<vec3> coeffs_u, coeffs_v, coeffs;
            compute_partial_derivative(this->coefficients, U, coeffs_u);
            compute_partial_derivative(this->coefficients, V, coeffs_v);
            cross(coeffs_u, coeffs_v, coeffs);
            result.set(coeffs);
        }

        Monomial2D &operator=(const Monomial2D &other) {
            if (this == &other) return *this; // Handle self-assignment
            n = other.n, m = other.m;
            coefficients = other.coefficients;
            return *this;
        }
    };

}
#endif //MONOMIAL_H
