//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <cmath>
#include <stdexcept>
#include <array>
#include <algorithm>
#include "cmmcore/binom.h"
namespace cmmcore{


class Matrix {
private:
    std::vector<double> data;
    size_t rows, cols;

public:
    Matrix(size_t r, size_t c) : rows(r), cols(c), data(r * c, 0.0) {}

    double& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }

    const double& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::runtime_error("Matrix dimensions mismatch for multiplication");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    Matrix inverse() const {
        if (rows != cols) {
            throw std::runtime_error("Only square matrices can be inverted");
        }
        Matrix augmented(rows, 2 * cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                augmented(i, j) = (*this)(i, j);
                augmented(i, j + cols) = (i == j) ? 1.0 : 0.0;
            }
        }

        for (size_t i = 0; i < rows; ++i) {
            size_t pivot = i;
            for (size_t j = i + 1; j < rows; ++j) {
                if (std::abs(augmented(j, i)) > std::abs(augmented(pivot, i))) {
                    pivot = j;
                }
            }
            if (pivot != i) {
                for (size_t j = 0; j < 2 * cols; ++j) {
                    std::swap(augmented(i, j), augmented(pivot, j));
                }
            }

            double pivotValue = augmented(i, i);
            if (std::abs(pivotValue) < 1e-10) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }

            for (size_t j = 0; j < 2 * cols; ++j) {
                augmented(i, j) /= pivotValue;
            }

            for (size_t k = 0; k < rows; ++k) {
                if (k != i) {
                    double factor = augmented(k, i);
                    for (size_t j = 0; j < 2 * cols; ++j) {
                        augmented(k, j) -= factor * augmented(i, j);
                    }
                }
            }
        }

        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = augmented(i, j + cols);
            }
        }
        return result;
    }
};

using Tensor3D = std::vector<Matrix>;

}


#endif //MATRIX_H
