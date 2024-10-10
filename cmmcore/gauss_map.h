//
// Created by Andrew Astakhov on 08.10.24.
//

#ifndef GAUSS_MAP_H
#define GAUSS_MAP_H

#ifdef CYTHON_ABI
#include "vec.h"
#include "convexhull.h"
#include "polygon.h"
#include "monomial.h"
#include "nurbs.h"
#else
#include "cmmcore/vec.h"
#include "cmmcore/convexhull.h"
#include "cmmcore/polygon.h"
#include "cmmcore/monomial.h"
#include "cmmcore/nurbs.h"
#endif
namespace cmmcore{

class GaussMap {


    public:
    NURBSSurface surf;
    std::vector<vec2> hull;

    GaussMap()=default;
    GaussMap( NURBSSurface& _surf, std::vector<vec2> &_hull):surf(std::move(_surf)),hull(std::move(_hull)) {

    };
    GaussMap(const NURBSSurface& srf) {
        build(srf);

      }
    void buildv2(const NURBSSurface& srf) {
      NURBSSurface a,b;
      srf.getDerivativeSurface(0,a);
      srf.getDerivativeSurface(1,b);
      surfaceCrossProduct(a,b,surf);

      std::vector<vec2> pts(surf._size[0]*surf._size[1]);
      auto flatcpts=surf.control_points_flat3d();
      for (int i = 0; i < (surf._size[0]*surf._size[1]); i++) {
        auto& pt=flatcpts[i];
        pt.unitize();
        cartesian_to_spherical(pt, pts[i]);
      }

      hull=convex_hull2d(pts);

    }
  void build(const NURBSSurface & srf) {

    //Tensor3D control_points(srf._size[0], Matrix(srf._size[1], 3));
      Monomial2D normal;
      Monomial2D mono=Monomial2D(surf);
      mono.computeNormal(normal);
      mono.to_bezier(surf);
      hull=std::vector<vec2>();
      std::vector<vec2> pts{};
      pts.resize(surf._size[0]*surf._size[1]);
      auto flatcpts=surf.control_points_flat3d();
      for (int i = 0; i < (surf._size[0]*surf._size[1]); i++) {
        auto& pt=flatcpts[i];
        pt.unitize();
        cartesian_to_spherical(pt, pts[i]);
      }

      hull=convex_hull2d(pts);


  }
  bool separabilityTest(const GaussMap& other) {

    switch (polygonRelationship2D(hull,other.hull)) {
      case PolygonRelationship2D::INTERSECT:
        return false;
      case PolygonRelationship2D::TOUCH:
        return true;
      case PolygonRelationship2D::DISTANCE:
        return true;
      default:
        throw std::runtime_error("GaussMap::separabilityTest: Unknown polygonRelationship");
    }


  }
    void subdivide(GaussMap& g11,GaussMap& g12,GaussMap& g21,GaussMap& g22)  {

      auto umid=0.5*(surf._interval[0][1]+surf._interval[0][0]);
      auto vmid=0.5*(surf._interval[1][1]+surf._interval[1][0]);
      auto [s1,s2]=surf.split_surface_u(       umid);
      auto [s11,s12]=s1.split_surface_v( vmid);
      auto [s21,s22]=s2.split_surface_v( vmid);

      std::vector<vec2> pts{};


      cartesian_to_spherical(s11.control_points_flat3d(),pts);

      g11.surf=std::move(s11);
      g11.hull=std::move(convex_hull2d(pts));

      pts.clear();

      cartesian_to_spherical(s12.control_points_flat3d(),pts);

      g12.surf=std::move(s12);
      g12.hull=std::move(convex_hull2d(pts));
      pts.clear();

      cartesian_to_spherical(s21.control_points_flat3d(),pts);

      g21.surf=std::move(s21);
      g21.hull=std::move(convex_hull2d(pts));

      pts.clear();

      cartesian_to_spherical(s22.control_points_flat3d(),pts);


      g22.surf=std::move(s22);
      g22.hull=std::move(convex_hull2d(pts));


    }



};

}
#endif //GAUSS_MAP_H
