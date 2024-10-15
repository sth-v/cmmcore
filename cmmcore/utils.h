/*
Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef UTILS_H
#define UTILS_H
#include <cmath>
#include <chrono>
#include <random>
#include <limits>
#include "cmmcore/vec.h"
namespace cmmcore {
#define EPS_DELTA 1e-12

    template <typename T>
    void vectorColumn(const std::vector<std::vector<T>>& matrix, const  size_t columnIndex,std::vector<T>& column) {
        column.resize(matrix[0].size());
        size_t i=0;
        for (const auto& row : matrix) {

            if (columnIndex < row.size()) {
                column[i]=row[columnIndex];
            } else {
                throw std::out_of_range("Column index is out of range");
            }
            i++;
        }

    }

/**
 * @brief Calculates an epsilon value for a given number.
 *
 * @param x The input number.
 * @return The calculated epsilon value.
 */

inline double calc_epsilon(const double x) noexcept {
    const double relative_epsilon = std::numeric_limits<double>::epsilon() * std::fabs(x);
    constexpr double absolute_epsilon = std::numeric_limits<double>::denorm_min();

    if ( std::fabs(x) < EPS_DELTA) {
        return absolute_epsilon;
    }
    return relative_epsilon;
}
    inline double random_double(double const min, double const max) {static std::random_device rd;static std::mt19937 gen(rd());std::uniform_real_distribution<> dis(min, max);return static_cast<double>(dis(gen));
}

    inline vec3 random_vec3(double const min, double const max) {
    return{random_double(min, max), random_double(min, max),random_double(min, max)};
}

    inline void random_vectors3(double const min, double const max, size_t const count, std::vector<vec3>& vecs) {
    vecs.resize(count);
    for (size_t i = 0; i < count; ++i) {
        vecs[i] = random_vec3(min, max);
    }
}
    class Timer {
    std::chrono::high_resolution_clock::time_point _start  ;
    std::chrono::high_resolution_clock::time_point  _end;
        public:
        int print_format=0;

        long long duration;
        Timer(int fmt=0) :_start(std::chrono::high_resolution_clock::now()),
        _end(_start),print_format(fmt),
        duration(0) {};
        void start() {
            this->_start=std::chrono::high_resolution_clock::now();

        }
        void stop() {
            this->_end=std::chrono::high_resolution_clock::now();
            this->duration=std::chrono::duration_cast<std::chrono::nanoseconds>(_end-_start).count();

        }
        void print(const std::string& text="") {
            if (print_format==0) {
                std::cout <<text<< duration << " ns." <<std::endl;
            }
            if (print_format==1) {

                    std::cout <<text<< duration*1e-6 << " ms." <<std::endl;
            }    if (print_format==2) {
                std::cout <<text<< duration*1e-9 << " secs." <<std::endl;

            }
        }

    };
}
#endif //UTILS_H
