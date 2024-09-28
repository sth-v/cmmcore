//
// Created by Andrew Astakhov on 19.09.24.
//
#include "cmmcore/nurbs.h"

#include <iostream>
using namespace cmmcore;
int main() {
  std::vector<double> U = {0, 0, 0, 0.5, 1, 1, 1};
  int n = 4;
  int p = 2;
  double u = 0.3;
  bool is_periodic = false;

  int span = find_span(n, p, u, U, is_periodic);
  std::cout << "The span index is: " << span << std::endl;  // Output: The span index is: 2

  return 0;
}

//int main() {
  // Example usage of the NURBS utilities

  // Assume we have n+1 control points in u-direction and m+1 in v-direction
  //int n = /* number of control points in u direction */ - 1;
  //int m = /* number of control points in v direction */ - 1;
  //int p = /* degree in u direction */;
  //int q = /* degree in v direction */;
  //
  //// Knot vectors
  //const double* U = /* pointer to knot vector in u direction */;
  //const double* V = /* pointer to knot vector in v direction */;
  //
  //// Control points (flattened array of Vec4)
  //const cmmcore::nurbs::Vec4* Pw = /* pointer to control points array */;
  //
  //// Parameter values
  //double u = 0.4/* parameter value in u */;
  //double v = 0.7/* parameter value in v */;

  // Initialize memory pool with sufficient size
  //cmmcore::nurbs::MemoryPool<double> pool(1024 * 1024); // Adjust size as needed

  // Evaluate point on NURBS surface
  //cmmcore::nurbs::Vec3 point;
  //cmmcore::nurbs::surface_point(n, p, U, m, q, V, Pw, u, v, false, false, point, pool);

  // Output the point
  //std::cout << "Point on NURBS surface: (" << point.x << ", " << point.y << ", " << point.z << ")\n";

//  return 0;
//}