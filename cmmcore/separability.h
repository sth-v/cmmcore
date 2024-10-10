//
// Created by Andrew Astakhov on 10.10.24.
//

#ifndef SEPARABILITY_H
#define SEPARABILITY_H
#include "cmmcore/vec.h"
namespace cmmcore {
/*
 * If the two objects are already known to intersect at a point and one
 * wishes to determine if they intersect in no other points then a variation
 * of 3-D separability test must be employed. Consider, for example
 * the situation depicted in Figure 4.4. A curve and a surface are known
 * to intersect at the point C. One wishes to know if they intersect at any
 * other point. The surface is split into four sub-patches and the curve into
 * two sub-curves. Every point except C is mapped to a point on the unit sphere
 * centered at C by a simple central projection. If a sub-patch and a sub-curve
 * intersect then their projections onto the sphere must intersect. If one
 * projects every point of the sub-curve onto the sphere and every point of the
 * sub-patch onto the sphere and these sets do not intersect, then the sub-curve
 * and sub-patch cannot intersect except at C. Practically, this is accomplished
 * as follows. If every control vertex of the Bézier or B-spline surface (curve)
 * except the control vertex corresponding to the surface corner (curve end) at C
 * is projected onto the sphere, the convex hull of the resulting point set contains
 * the projection of the surface (curve) onto the sphere.
 *
 * Если известно, что два объекта пересекаются в одной точке, и нужно определить, пересекаются ли они в других точках,
 * то необходимо использовать разновидность теста на трехмерную делимость. Рассмотрим, например, ситуацию, изображенную
 * на рисунке 4.4. Известно, что кривая и поверхность пересекаются в точке C. Требуется определить, пересекаются ли они
 * в любой другой точке. Поверхность разбивается на четыре подпатча, а кривая - на две подкривые. Каждая точка, кроме C,
 * отображается на точку единичной сферы с центром в C с помощью простой центральной проекции. Если подпатч и подкривая
 * пересекаются, то их проекции на сферу должны пересекаться. Если спроецировать каждую точку подкривой на сферу и
 * каждую точку подпатча на сферу и эти множества не пересекаются, то подкривая и подпатч не могут пересекаться, кроме
 * как в точке C. Практически это достигается следующим образом. Если каждую контрольную вершину поверхности (кривой)
 * Безье или B-сплайна, кроме контрольной вершины, соответствующей углу поверхности (концу кривой) в точке C,
 * спроецировать на сферу, то выпуклая оболочка полученного множества точек содержит проекцию поверхности (кривой) на
 * сферу.
 * */

class SphericalBBox{
  public:
    vec3 start;
    vec3 end;

    SphericalBBox()=default;

    SphericalBBox(const std::vector<vec3>& pts, const vec3& center=0):start(0),end(0){

      std::vector<vec3> sphere_points1;
      project_to_sphere(pts, center, sphere_points1);
      for (int axis = 0; axis < 3; ++axis) {
        std::vector<vec2> plane_coords1;
        plane_coords1.reserve(sphere_points1.size());

        for (const auto& pt : sphere_points1) {
          if (axis == 0) {
            plane_coords1.emplace_back(pt.y, pt.z);
          } else if (axis == 1) {
            plane_coords1.emplace_back(pt.z, pt.x);
          } else {
            plane_coords1.emplace_back(pt.x, pt.y);
          }
        }
        find_smallest_wedge(plane_coords1, axis);
      }


    }
  bool separable(const SphericalBBox& other) const {
      for (int axis = 0; axis < 3; ++axis) {

        if ((end[axis]< other.start[axis] && other.end[axis] > start[axis]) ||     (other.end[axis] < start[axis] && end[axis] > other.start[axis])){
          return true;
        }
      }
      return false;
    };
  bool intersects(const SphericalBBox& other) const {
    return !separable(other);
  }
  private:
  void project_to_sphere(const std::vector<vec3>& points,const vec3& center, std::vector<vec3>& result) {
      result.resize(points.size());
      for ( int i=0;i<points.size();++i ) {
        const auto& point=points[i];
        point.sub(center,  result[i]);

        result[i].unitize();


      }

    }

  void find_smallest_wedge(const std::vector<vec2>& points_2d,size_t dim) {
      std::vector<double> angles(points_2d.size());
      for (int i = 0; i < points_2d.size(); ++i) {
        angles[i]=std::atan2(points_2d[i].y, points_2d[i].x);
      }
      std::sort(angles.begin(), angles.end());

      std::vector<double> angle_diffs(angles.size());
      std::adjacent_difference(angles.begin(), angles.end(), angle_diffs.begin());
      angle_diffs[0] = angles[0] + 2 * M_PI - angles.back(); // Handle circular difference

      auto max_gap_it = std::max_element(angle_diffs.begin(), angle_diffs.end());
      size_t max_gap_idx = std::distance(angle_diffs.begin(), max_gap_it);

      start[dim]= angles[max_gap_idx];
      end[dim] = angles[(max_gap_idx + 1) % angles.size()];



    }


};


}
#endif //SEPARABILITY_H
