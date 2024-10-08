//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef MONOMIAL_H
#define MONOMIAL_H

#ifdef CYTHON_ABI
#include "matrix.h"
#include "binom.h"
#else
#include "cmmcore/matrix.h"
#include "cmmcore/binom.h"
#endif

namespace cmmcore {
inline Matrix bpmat(int n) {
    Matrix result(n + 1, n + 1);
    for (int j = 0; j <= n; ++j) {
        for (int k = 0; k <= j; ++k) {
            result(j, k) = std::pow(-1, j - k) * binomial_coefficient(n, j) * binomial_coefficient(j, k);
        }
    }
    return result;
}

inline Tensor3D bezier_to_monomial(const Tensor3D& control_points) {
    size_t n = control_points.size();
    size_t m = control_points[0].getRows();
    size_t l = control_points[0].getCols();

    Matrix Mu = bpmat(n - 1);
    Matrix Mv = bpmat(m - 1);

    Tensor3D monomial_coeffs(l, Matrix(n, m));
    for (size_t d = 0; d < l; ++d) {
        Matrix cp_slice(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                cp_slice(i, j) = control_points[i](j, d);
            }
        }
        monomial_coeffs[d] = Mu * cp_slice * Mv.transpose();
    }

    return monomial_coeffs;
}

inline Tensor3D monomial_to_bezier(const Tensor3D& monomial_coeffs) {
    size_t n = monomial_coeffs[0].getRows();
    size_t m = monomial_coeffs[0].getCols();
    size_t l = monomial_coeffs.size();

    Matrix Mu_inv = bpmat(n - 1).inverse();
    Matrix Mv_inv = bpmat(m - 1).inverse();

    Tensor3D control_points(n, Matrix(m, l));
    for (size_t d = 0; d < l; ++d) {
        Matrix bez_slice = Mu_inv * monomial_coeffs[d] * Mv_inv.transpose();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                control_points[i](j, d) = bez_slice(i, j);
            }
        }
    }

    return control_points;
}

inline  Tensor3D cross_product_monomial(const Tensor3D& a_coeffs, const Tensor3D& b_coeffs) {
    size_t n1 = a_coeffs[0].getRows(), m1 = a_coeffs[0].getCols();
    size_t n2 = b_coeffs[0].getRows(), m2 = b_coeffs[0].getCols();
    size_t n = n1 + n2 - 1, m = m1 + m2 - 1;

    Tensor3D result(3, Matrix(n, m));

    for (size_t i1 = 0; i1 < n1; ++i1) {
        for (size_t j1 = 0; j1 < m1; ++j1) {
            for (size_t i2 = 0; i2 < n2; ++i2) {
                for (size_t j2 = 0; j2 < m2; ++j2) {
                    size_t i = i1 + i2, j = j1 + j2;
                    result[0](i, j) += a_coeffs[1](i1, j1) * b_coeffs[2](i2, j2) - a_coeffs[2](i1, j1) * b_coeffs[1](i2, j2);
                    result[1](i, j) += a_coeffs[2](i1, j1) * b_coeffs[0](i2, j2) - a_coeffs[0](i1, j1) * b_coeffs[2](i2, j2);
                    result[2](i, j) += a_coeffs[0](i1, j1) * b_coeffs[1](i2, j2) - a_coeffs[1](i1, j1) * b_coeffs[0](i2, j2);
                }
            }
        }
    }

    return result;
}

inline std::pair<Tensor3D, Tensor3D> monomial_partial_derivatives(const Tensor3D& coeffs) {
    size_t n = coeffs[0].getRows(), m = coeffs[0].getCols(), dim = coeffs.size();

    Tensor3D du_coeffs(dim, Matrix(n - 1, m));
    Tensor3D dv_coeffs(dim, Matrix(n, m - 1));

    for (size_t d = 0; d < dim; ++d) {
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                du_coeffs[d](i - 1, j) = i * coeffs[d](i, j);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 1; j < m; ++j) {
                dv_coeffs[d](i, j - 1) = j * coeffs[d](i, j);
            }
        }
    }

    return {du_coeffs, dv_coeffs};
}

inline Tensor3D normal_vector_monomial(const Tensor3D& coeffs) {
    auto [du_coeffs, dv_coeffs] = monomial_partial_derivatives(coeffs);
    return cross_product_monomial(du_coeffs, dv_coeffs);
}
}
#endif //MONOMIAL_H
