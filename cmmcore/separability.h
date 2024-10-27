//
// Created by Andrew Astakhov on 10.10.24.
//

#ifndef SEPARABILITY_H
#define SEPARABILITY_H
#include "cmmcore/vec.h"
#include "cmmcore/aabb.h"
#include "cmmcore/polygon.h"
#include "cmmcore/convexhull.h"
#include <algorithm>
#include "cmmcore/gjk.h"
#include "cmmcore/numeric_utils.h"
#include "cmmcore/utils.h"
#ifdef CMMCORE_DEBUG

#endif
namespace cmmcore
{

  inline bool AABBTest(const std::vector<vec3>& pts1,const std::vector<vec3>& pts2) {
    const AABB bb1(pts1);
    const AABB bb2(pts2);
    return !bb1.intersects(bb2);
  }
  inline void polygonProject(const polygon2 &polygon, const vec2 &direction, vec2 &minmax) {
    minmax.x = std::numeric_limits<double>::max();
    minmax.y = std::numeric_limits<double>::lowest();
    for (auto i: polygon) {
      auto res = i.projection(direction);
      minmax.x = std::min(minmax.x, res);
      minmax.y = std::max(minmax.y, res);
    }
  };
  ;

  inline std::pair<vec2, vec2> swapIntervals(vec2 &mm1, vec2 &mm2) {
    if (mm1.y > mm2.y) {
      return {mm2, mm1};
    } else { return {mm1, mm2}; }
  };

  inline bool SAT2D(const polygon2 &a, const polygon2 &b) {

    UnorderedNormal2Set normals;
    if (a.empty()|| b.empty()){
      throw std::runtime_error("SAT2D: empty polygon");
    }
    if (b.size()<2 && a.size()<2)
    {
      return b!=a;
    }
    if (a.size()<2)
    {
      for (size_t i = 0; i <b.size(); i++)
      {
        size_t j = (i + 1) % b.size();
        vec2 n = {b[j].y - b[i].y, b[i].x - b[j].x};

        if (contains(normals, n)) {
          continue;
        }
        vec2 minmax1, minmax2;

        polygonProject(a, n, minmax1);
        polygonProject(b, n, minmax2);

        if (auto [m1,m2] = swapIntervals(minmax1, minmax2); m1.y < m2.x) {
          return true;
        };


        normals.insert(n);
      }
    } else
    {
      for (size_t i = 0; i < a.size(); i++)
      {
        size_t j = (i + 1) % a.size();
        vec2 n = {a[j].y - a[i].y, a[i].x - a[j].x};

        if (contains(normals, n)) {
          continue;
        }
        vec2 minmax1, minmax2;
        polygonProject(a, n, minmax1);
        polygonProject(b, n, minmax2);

        if (auto [m1,m2] = swapIntervals(minmax1, minmax2); m1.y < m2.x) {
          return true;
        };


        normals.insert(n);
      }
    }
    return false;
  }
  /*
    class  SAT3DHandler {
      const std::vector<vec3>* pts1;
      const std::vector<vec3>*  pts2;
      const double eps;
      std::vector<vec3> hull1,hull2;
      ConvexHullResult chtype1;
      ConvexHullResult chtype2;
      // true - separable, false - not separable

      static bool handleEquality( std::vector<vec3>& h1, std::vector<vec3>& h2) {
        if( &h1 == &h2) {
          return false;
        }

        for (size_t i=0;i<h1.size();i++) {
          if (h1[i]==h2[i]) {
            return false;
          }
        }
        return false;


      }
      static bool handlePointPoint( std::vector<vec3>& h1, std::vector<vec3>& h2) {
        return h1[0]==h2[0];

      }

      static bool handlePointLine( std::vector<vec3>& h1, std::vector<vec3>& h2) {
        const auto d2=(h2[1]-h2[0]);
        const auto d1=(h1[0]-h2[0]);
        const auto d1p=d1.project(d2);
        return (d1p-d1).length()<=std::numeric_limits<double>::epsilon();


      }
       static bool handlePointFlat( std::vector<vec3>& h1, std::vector<vec3>& h2) {
        plane p={{h2[0],h2[1],h2[2]}};
        double val=p.signedDistance(h1[0]);

        if (std::abs(val)<=std::numeric_limits<double>::epsilon()) {
          polygon2 h22d(h2.size());
          vec2 pt;

          if (p.normal.dot({0,0,1})!=0) {
            pt.set(h1[0].x,h1[0].y);
            for (size_t i=0;i<h2.size();i++) {
              h22d[i].set(h2[i].x,h2[i].y);


            }
          } else if (p.normal.dot({0,1,0})!=0) {
            pt.set(h1[0].x,h1[0].z);
            for (size_t i=0;i<h2.size();i++) {
              h22d[i].set(h2[i].x,h2[i].z);

            }
          }else {
            pt.set(h1[0].z,h1[0].y);
            for (size_t i=0;i<h2.size();i++) {
              h22d[i].set(h2[i].z,h2[i].y);

            }
          }
          return isvec2Inside2d(h22d,pt);



          }
          return false;






      }
       static bool handlePointVolume( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }
       static bool handleLineLine( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }

       static bool handleLineFlat( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }
       static bool handleLineVolume( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }

       static bool handleFlatFlat( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }
      static  bool handleFlatVolume( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }
       static bool handleVolumVolume( std::vector<vec3>& h1, std::vector<vec3>& h2) {


      }
    public:


       SAT3DHandler(const std::vector<vec3>& _pts1,const std::vector<vec3>& _pts2, const double _eps=std::numeric_limits<double>::epsilon()):pts1(&_pts1),pts2(&_pts2), eps(_eps){

      }
      bool solve() {
        chtype1=convex_hull(*pts1,hull1);
        chtype2=convex_hull(*pts2,hull2);
        return false;
      }


    };
    */
  inline bool SAT3D(const std::vector<vec3>& pts1,const std::vector<vec3>& pts2, const double eps=std::numeric_limits<double>::epsilon()) {
    std::vector<vec3> hull1,hull2;
#ifdef CMMCORE_DEBUG
    auto t=Timer(1);
    t.start();
#endif

    convex_hull3d(pts1,hull1);
#ifdef CMMCORE_DEBUG
    t.stop();
    t.print("convex_hull3d at: ");
    t.start();
#endif


    convex_hull3d(pts2,hull2);
#ifdef CMMCORE_DEBUG
    t.stop();
    t.print("convex_hull3d at: ");
#endif

    std::vector<vec3>  simplex;
    vec3 closestPointToOrigin;
#ifdef CMMCORE_DEBUG
    t.start();
#endif

    GJK(hull1,hull2, simplex,closestPointToOrigin,eps);
#ifdef CMMCORE_DEBUG
    auto ff=closestPointToOrigin.length();
    t.stop();
    printf("%f ",ff );
    ff>eps? printf("no collision"): printf("collision detected");
    t.print("\nGJK at: ");


#else
    auto ff=closestPointToOrigin.length();
#endif


    return ff>eps;
  }







  /*
   * SPHERICAL SEPARABILITY
   *
   *
   */



  inline bool SphericalAABBTest( std::vector<vec3>& pts1, std::vector<vec3>& pts2) {
    const SphericalAABB bb1(pts1);
    const SphericalAABB bb2(pts2);
    return bb1.separable(bb2);
  }

  using Vec3Set=std::unordered_set<vec3, Vec3Hash,Vec3Equal> ;


  inline bool contains( const Vec3Set& s, const vec3& val) {
    return s.find(val) != s.end();
  }



  inline bool SphericalCentralProjectionTest(std::vector<vec3>& pts1, std::vector<vec3>& pts2) {

    std::vector<vec2> poly1(pts1.size());
    std::vector<vec2> poly2(pts2.size());

    vec3 O(0,0,1);
    for (size_t i=0;i<pts1.size();i++)
    {

      centerProjectionToXY(O, pts1[i].unit(),poly1[i]);
    }
    for (size_t j=0;j<pts2.size();j++)
    {
      centerProjectionToXY(O, pts2[j].unit(),poly2[j]);

    }

#ifdef CMMCORE_DEBUG
    printf("\n\n[");
    for (auto& p:poly1)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("],[");
    for (auto& p:poly2)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("]\n");
#endif

    auto h1=convex_hull2d(poly1);
     auto h2=convex_hull2d(poly2);
    if (h1.empty())
    {
      h1.push_back(poly1[0]);
    }
    if (h2.empty())
    {
      h2.push_back(poly2[0]);
    }
#ifdef CMMCORE_DEBUG
    printf("[");
    for (auto& p:h1)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("],[");
    for (auto& p:h2)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("]\n");
#endif

    return SAT2D(h1,h2);

    }
inline bool SphericalSeparabilityTest( std::vector<vec3>& pts1, std::vector<vec3>& pts2) {

    //printf(("SEPTEST\n\n["+format_vec3vec(pts1)+"],["+format_vec3vec(pts2)+"]").c_str());
   // if (SphericalAABBTest(pts1,pts2))
    //{

   //   printf("BBB");
   //   return true;
    //}
    return SphericalCentralProjectionTest(pts1,pts2);

  }

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

  inline bool SphericalCentralProjectionTest2(vec3 origin, std::vector<vec3>& pts1, std::vector<vec3>& pts2) {

    std::vector<vec2> poly1(pts1.size());
    std::vector<vec2> poly2(pts2.size());


    for (size_t i=0;i<pts1.size();i++)
    {

      centerProjectionToXY(origin, pts1[i].unit(),poly1[i]);
    }
    for (size_t j=0;j<pts2.size();j++)
    {
      centerProjectionToXY(origin, pts2[j].unit(),poly2[j]);

    }

#ifdef CMMCORE_DEBUG
    printf("\n\n[");
    for (auto& p:poly1)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("],[");
    for (auto& p:poly2)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("]\n");
#endif

    auto h1=convex_hull2d(poly1);
    auto h2=convex_hull2d(poly2);
    if (h1.empty())
    {
      h1.push_back(poly1[0]);
    }
    if (h2.empty())
    {
      h2.push_back(poly2[0]);
    }
#ifdef CMMCORE_DEBUG
    printf("[");
    for (auto& p:h1)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("],[");
    for (auto& p:h2)
    {
      printf("[%f,%f],", p.x,p.y);
    }
    printf("]\n");
#endif

    return SAT2D(h1,h2);

  }
  }

#endif //SEPARABILITY_H
