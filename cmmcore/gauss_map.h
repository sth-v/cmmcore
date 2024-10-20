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
#include "separability.h"

namespace cmmcore
{
  inline void inverseControlPoints(  const std::vector<vec3>& pts, std::vector<vec3>&pts_inv)
  {
    pts_inv.resize(pts.size());



    for (size_t i = 0; i < pts.size(); ++i)
    {
      pts_inv[i].set(-pts[i].x, -pts[i].y, -pts[i].z);

    }
  }
  class GaussMap
  {


  public:
    std::vector<vec3> pts;
    std::vector<vec3> pts_inv;

    void solveInvPts()
    {
      pts_inv.clear();
      pts_inv.reserve(pts.size());
      for (auto& p : pts)
      {
        pts_inv.emplace_back(-p.x, -p.y, -p.z);

      }
    }

    NURBSSurface surf{};


    GaussMap() = default;

    GaussMap( NURBSSurface& srf)
    {
      build(srf);
    }


    void build( NURBSSurface& srf)
    {
      //Tensor3D control_points(srf._size[0], Matrix(srf._size[1], 3));

      auto mono = Monomial2D(srf);
      Monomial2D normal;
      mono.computeNormal(normal);
      normal.to_bezier(surf);
      surf.control_points_flat3d(pts);


      solveInvPts();

      //cartesian_to_spherical(flatcpts[i], pts[i]);
    }

    bool separabilityTest(GaussMap& other)
    {

      SphericalSeparabilityTest(pts, other.pts);
      return (SphericalSeparabilityTest(pts, other.pts) &&
        SphericalSeparabilityTest(pts, other.pts_inv) &&
        SphericalSeparabilityTest(pts_inv, other.pts_inv) &&
        SphericalSeparabilityTest(pts_inv, other.pts));
    }


    void subdivide(GaussMap& g11, GaussMap& g12, GaussMap& g21, GaussMap& g22)
    {
      const auto pts1 = surf.control_points_flat3d();
      auto umid = 0.5 * (surf._interval[0][1] + surf._interval[0][0]);
      auto vmid = 0.5 * (surf._interval[1][1] + surf._interval[1][0]);
      surf.subdivide(g11.surf,g12.surf,g21.surf,g22.surf, umid, vmid);

      //auto [s1,s2] = surf.split_surface_u(umid);
      //auto [s11,s12] = s1.split_surface_v(vmid);
      //auto [s21,s22] = s2.split_surface_v(vmid);

       g11.surf.control_points_flat3d( g11.pts);
      g11.solveInvPts();


      g12.surf.control_points_flat3d( g12.pts);
      g12.solveInvPts();


      g21.surf.control_points_flat3d( g21.pts);
      g21.solveInvPts();


      g22.surf.control_points_flat3d( g22.pts);
      g22.solveInvPts();


    }
  };
  inline NURBSSurface gaussMap(const NURBSSurface& srf)
  {
    auto mono = Monomial2D(srf);
    Monomial2D normal;
    mono.computeNormal(normal);
    NURBSSurface surf{};
    normal.to_bezier(surf);
    return surf;
  }
  inline bool gaussMapSeparability(const NURBSSurface& srf1,const NURBSSurface& srf2)
  {
    std::vector<vec3>pts1=srf1.control_points_flat3d();;
    std::vector<vec3>pts2=srf2.control_points_flat3d();;
    std::vector<vec3>pts1_inv;
    std::vector<vec3>pts2_inv;
    inverseControlPoints(pts1,pts1_inv);
    inverseControlPoints(pts2,pts2_inv);

    return (SphericalSeparabilityTest(pts1,pts2)&&
    SphericalSeparabilityTest(pts1,pts2_inv)&&
    SphericalSeparabilityTest(pts1_inv,pts2_inv)&&
    SphericalSeparabilityTest(pts1_inv,pts2));


  }
}
#endif //GAUSS_MAP_H
