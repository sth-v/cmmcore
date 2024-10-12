//
// Created by Andrew Astakhov on 11.10.24.
//

#ifndef PLANE_H
#define PLANE_H
#include "cmmcore/vec.h"
#include <tuple>
namespace cmmcore {
    inline void plane_eq_from_pts(const vec3 &pt0, const vec3 &pt1, const vec3 &pt2, vec4& equation) {
        const float x3 = pt0.y * pt1.z;
        const float x4 = pt1.y * pt2.z;
        const float x5 = pt2.y * pt0.z;
        equation.x = x3 + x4 + x5 - pt0.y * pt2.z - pt1.y * pt0.z - pt2.y * pt1.z;
        equation.y = -pt0.x * pt1.z + pt0.x * pt2.z + pt1.x * pt0.z - pt1.x * pt2.z - pt2.x * pt0.z + pt2.x * pt1.z;
        equation.z = pt0.x * pt1.y - pt0.x * pt2.y - pt1.x * pt0.y + pt1.x * pt2.y + pt2.x * pt0.y - pt2.x * pt1.y;
        equation.w = -pt0.x * x4 + pt0.x * pt2.y * pt1.z - pt1.x * x5 + pt1.x * pt0.y * pt2.z - pt2.x * x3 + pt2.x * pt1.y * pt0.z;

    }
    inline double planeSignedDistance(const vec3 &pt, const vec4 &plane) {
        return plane.x*pt.x + plane.y*pt.y + plane.z*pt.z + plane.w;
    };
    class plane{
        double denominator=1;
    public:

        vec4 equation{0,0,1,0};
        vec3 normal{0,0,1};
        vec3 origin{0,0,0};
        plane()=default;
        plane(double a,double b,double c,double d=0) : equation(a,b,c,d) {
            computeDenominator();

            normal = vec3(equation.x/denominator,equation.y/denominator,equation.z/denominator);
            origin.set(normal*equation.w);
        }
        plane(const vec3& _normal) : denominator(_normal.length()), equation(_normal.x,_normal.y,_normal.z, 0),normal(_normal/denominator)  {
        }
        plane(const vec3& _origin,const vec3& _normal) : denominator(_normal.length()), equation(_normal.x,_normal.y,_normal.z, -(_normal.x*_origin.x,_normal.y*_origin.y,_normal.z*_origin.z)), normal(_normal/denominator),origin(_origin)  {
        }

        plane(const std::tuple<vec3,vec3,vec3>& _pts)  {
            plane_eq_from_pts(std::get<0>(_pts),std::get<1>(_pts),std::get<2>(_pts),equation);
            computeDenominator();
            normal = vec3(equation.x/denominator,equation.y/denominator,equation.z/denominator);
            origin.set(normal*equation.w);
        };
        void computeDenominator() {
            denominator = std::sqrt(equation.x*equation.x+equation.y*equation.y+equation.z*equation.z);

        }
        double signedDistance(const vec3& pt) const {
            return (pt.x*equation.x+pt.y*equation.y+pt.z*equation.z+equation.w)/denominator;
        };
        void signedDistance(const std::vector<vec3>& pts, std::vector<double>& result) const {
            result.resize(pts.size());
            for (size_t i = 0; i < pts.size(); ++i){
                result[i]=signedDistance(pts[i]);
            }
        }

    };
    inline bool isPtsCoplanar(const std::vector<vec3>& pts,const double eps=std::numeric_limits<double>::epsilon()) {
        if (pts.size()<=3) {
            return true;
        }
        plane p={{pts[0],pts[1],pts[2]}};


        for (size_t i = 4; i < pts.size(); i++) {


            if (std::abs(p.signedDistance(pts[i]))>=eps) {
                return false;

            }
        }
        return true;
    }

    class cplane: plane{
      public:

        vec3 xaxis{1,0,0};
        vec3 yaxis{0,1,0};

        cplane()=default;

        cplane(const vec3& _origin,const vec3& xaxis,const vec3& yaxis) : plane(_origin,xaxis.cross(yaxis)),xaxis(xaxis),yaxis(yaxis)  {

        }
        cplane(const std::tuple<vec3,vec3,vec3>& _pts): plane(_pts){
            xaxis=std::get<0>(_pts)-origin;
            xaxis.unitize();
            yaxis=normal.cross(xaxis);
            yaxis.unitize();
        };
        double signedDistance(const vec3& pt) const {return plane::signedDistance(pt-origin);};

        void evaluate(double u,double v, vec3& result) const {
            result.set(origin);
            result+=  (xaxis*u);
            result+=(yaxis*v);

        }
        void evaluateInv( const vec3& xyz, vec3& uvh) const {
            auto temp=xyz -origin;
            uvh.x=temp.dot(xaxis);
            uvh.y=temp.dot(yaxis);
            uvh.z=temp.dot(normal);
        }
        void evaluateInv( const std::vector<vec3>& xyz, std::vector<vec3>& uvh) const {
          uvh.resize(xyz.size());
          vec3 temp;
          for (size_t i = 0; i < xyz.size(); ++i){
            xyz[i].sub(origin,temp);
            uvh[i].x=temp.dot(xaxis);
            uvh[i].y=temp.dot(yaxis);
             uvh[i].z=temp.dot(normal);

          }


        }
        void planeLineIntersection(const vec3& p1,const vec3& p2, double& u, double& v) {
            // Implementation of the algorithm from https://www.geometrictools.com/Documentation/IntersectionLine3Plane3.pdf
            //
        }
        void planePlaneIntersection(const plane& other, vec3& p1, vec3& p2 ) {
            //...
        }
        void planePlanePlaneIntersection(const plane& other1,const plane& other2, vec3& pt ) {
            //...
        }
   };
};
#endif //PLANE_H
