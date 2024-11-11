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

namespace cmmcore{
    // Square matrix N x N template
    template <size_t N>
    using MatrixSq = std::array<std::array<double, N>, N>;

    using Matrix3x3=MatrixSq<3>;
    using Matrix4x4=MatrixSq<4>;

    // Matrix N x M size template
    template <size_t Rows, size_t Cols>
    using Matrix = std::array<std::array<double,Cols>, Rows>;

#include <array>
#include <stdexcept>

class Matrix2x2
{
private:
    std::array<std::array<double, 2>, 2> _m{};
public:
    Matrix2x2() = default;
    Matrix2x2( std::initializer_list<std::initializer_list<double>>arr){};
    Matrix2x2(std::array<double,2>& ab, std::array<double,2>& cd)
    {
        _m[0][0]=ab[0];
        _m[0][1]=ab[1];
        _m[1][0]=cd[0];
        _m[1][1]=cd[1];
    }
    Matrix2x2(const double a, const double b, const double c, const double d)
    {
        _m[0][0] = a;
        _m[0][1] = b;
        _m[1][0] = c;
        _m[1][1] = d;
    };
    Matrix2x2(const Matrix2x2& other) = default;

    static Matrix2x2 fromDiag(const std::array<double, 2>& diag)
    {
        Matrix2x2 mtrx{};
        mtrx._m[0][0] = diag[0];
        mtrx._m[1][1] = diag[1];
        return mtrx;
    }

    static Matrix2x2 eye()
    {
        Matrix2x2 mtrx{};
        mtrx._m[0][0] = 1.0;
        mtrx._m[1][1] = 1.0;
        return mtrx;
    }

    std::array<double, 2>& operator[](const size_t i)
    {
        return _m[i];
    }

    const std::array<double, 2>& operator[](const size_t i) const
    {
        return _m[i];
    }

    Matrix2x2 transpose() const {
        Matrix2x2 result;

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                result._m[j][i] = _m[i][j];
            }
        }
        return result;
    }

    // Compute the determinant of the matrix
    double det() const {
        return (_m[0][0] * _m[1][1]) - (_m[0][1] * _m[1][0]);
    }

    // Compute the inverse of the matrix
    Matrix2x2 inverse() const {
        double determinant = det();

        if (determinant == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        Matrix2x2 inv;
        inv._m[0][0] =  _m[1][1] / determinant;
        inv._m[0][1] = -_m[0][1] / determinant;
        inv._m[1][0] = -_m[1][0] / determinant;
        inv._m[1][1] =  _m[0][0] / determinant;

        return inv;
    }

    // Matrix multiplication
    Matrix2x2 operator*(const Matrix2x2& other) const {
        Matrix2x2 result;

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                result._m[i][j] = 0.0;
                for (size_t k = 0; k < 2; ++k) {
                    result._m[i][j] += _m[i][k] * other._m[k][j];
                }
            }
        }

        return result;
    }

    // Optional: Matrix-vector multiplication
    std::array<double, 2> operator*(const std::array<double, 2>& vec) const {
        std::array<double, 2> result = {0.0, 0.0};
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                result[i] += _m[i][j] * vec[j];
            }
        }
        return result;
    }
};


    // Solves the system A * x = b for x, where A is a 2x2 matrix and b is a 2D vector
    inline bool solve2x2(const Matrix2x2& A, const std::array<double, 2>& b, std::array<double, 2>& result) noexcept {
        const double detA = A.det();

        if (detA == 0.0) {
            return false;
        }

        // Create matrices for Cramer's Rule
        Matrix2x2 A1 = A;
        A1[0][0] = b[0];
        A1[1][0] = b[1];
        double detA1 = A1.det();

        Matrix2x2 A2 = A;
        A2[0][1] = b[0];
        A2[1][1] = b[1];
        const double detA2 = A2.det();

        // Compute the solution using Cramer's Rule
        result[0] = detA1 / detA;
        result[1]  = detA2 / detA;
        return true;


    }


    // Matrix N x M of dynamical size
    class MatrixD {
    public:
        std::vector<double> data;
        size_t rows=0, cols=0;

        MatrixD()=default;
        MatrixD(const size_t r, const size_t c) : rows(r), cols(c), data(r * c, 0.0) {}

        double& operator()(size_t i, size_t j) {
            return data[i * cols + j];
        }

        const double& operator()(size_t i, size_t j) const {
            return data[i * cols + j];
        }

        size_t  getRows() const { return rows; }
         size_t getCols() const { return cols; }

        MatrixD operator*(const MatrixD& other) const {
            if (cols != other.rows) {
                throw std::runtime_error("Matrix dimensions mismatch for multiplication");
            }
            MatrixD result(rows, other.cols);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < other.cols; ++j) {
                    for (size_t k = 0; k < cols; ++k) {
                        result(i, j) += (*this)(i, k) * other(k, j);
                    }
                }
            }
            return result;
        }

        MatrixD transpose() const {
            MatrixD result(cols, rows);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    result(j, i) = (*this)(i, j);
                }
            }
            return result;
        }

        MatrixD inverse() const {
            if (rows != cols) {
                throw std::runtime_error("Only square matrices can be inverted");
            }
            MatrixD augmented(rows, 2 * cols);
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

            MatrixD result(rows, cols);
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    result(i, j) = augmented(i, j + cols);
                }
            }
            return result;
        }
    };

}


#endif //MATRIX_H


