//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef CMMCORE_LP_H
#define CMMCORE_LP_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>

namespace cmmcore {
    namespace lp{


        typedef std::vector<double> Vector;

        // Constraint structure
        struct Constraint {
            Vector a; // Coefficients of the constraint
        };

        // Function to compute dot product
        double dot(const Vector& a, const Vector& b) {
            double result = 0.0;
            for (size_t i = 0; i < a.size(); i++)
                result += a[i] * b[i];
            return result;
        }

        // Function to subtract two vectors
        Vector subtract(const Vector& a, const Vector& b) {
            Vector result(a.size());
            for (size_t i = 0; i < a.size(); i++)
                result[i] = a[i] - b[i];
            return result;
        }

        // Function to add two vectors
        Vector add(const Vector& a, const Vector& b) {
            Vector result(a.size());
            for (size_t i = 0; i < a.size(); i++)
                result[i] = a[i] + b[i];
            return result;
        }

        // Function to multiply vector by scalar
        Vector scalarMultiply(const Vector& a, double s) {
            Vector result(a.size());
            for (size_t i = 0; i < a.size(); i++)
                result[i] = a[i] * s;
            return result;
        }

        // Linear Programming Solver class
        class LPSolver {
        public:
            // Solve the LP problem
            Vector solve(int dimension, const std::vector<Constraint>& constraints, const Vector& n, const Vector& d) {
                std::vector<Constraint> randomized_constraints = constraints;
                // Randomize constraints except for the first one (x_d >= 0)
                std::shuffle(randomized_constraints.begin() + 1, randomized_constraints.end(), rng);

                // Initialize provisional optimum
                Vector x(dimension + 1, 0.0);
                x[dimension] = 1.0; // Start with x_d = 1

                // Process constraints incrementally
                return solveRecursive(dimension, randomized_constraints, n, d, x, 0);
            }

        private:
            std::mt19937 rng{std::random_device{}()};

            // Recursive solver function
            Vector solveRecursive(int dimension, const std::vector<Constraint>& constraints, const Vector& n, const Vector& d, Vector x, int k) {
                if (k == constraints.size())
                    return x;

                Constraint a_k = constraints[k];

                // Check if x satisfies the constraint
                if (dot(a_k.a, x) >= 0) {
                    return solveRecursive(dimension, constraints, n, d, x, k + 1);
                } else {
                    // Need to find new provisional optimum
                    if (dimension == 0) {
                        // No feasible solution
                        return x;
                    }

                    // Build new problem with constraint a_k.a * x = 0
                    int new_dimension = dimension - 1;

                    // Build new constraints
                    std::vector<Constraint> new_constraints;
                    for (int i = 0; i < k; i++) {
                        Constraint c_i = constraints[i];

                        // Compute projection
                        double denom = dot(a_k.a, a_k.a);
                        if (denom == 0) continue; // Avoid division by zero

                        double factor = dot(c_i.a, a_k.a) / denom;
                        Constraint new_c;
                        new_c.a = subtract(c_i.a, scalarMultiply(a_k.a, factor));
                        new_constraints.push_back(new_c);
                    }

                    // Build new objective function
                    double denom = dot(a_k.a, a_k.a);
                    if (denom == 0) denom = 1e-8; // Avoid division by zero

                    double factor_n = dot(n, a_k.a) / denom;
                    Vector new_n = subtract(n, scalarMultiply(a_k.a, factor_n));

                    double factor_d = dot(d, a_k.a) / denom;
                    Vector new_d = subtract(d, scalarMultiply(a_k.a, factor_d));

                    // Recursive call
                    Vector x_new = solveRecursive(new_dimension, new_constraints, new_n, new_d, x, 0);

                    // Lift solution back to original dimension
                    double t_num = dot(a_k.a, x_new);
                    double t_denom = dot(a_k.a, x);
                    if (t_denom == 0) t_denom = 1e-8; // Avoid division by zero

                    double t = -t_num / t_denom;
                    x = add(x_new, scalarMultiply(x, t));
                    return solveRecursive(dimension, constraints, n, d, x, k + 1);
                }
            }
        };



    }
}
#endif //CMMCORE_LP_H
