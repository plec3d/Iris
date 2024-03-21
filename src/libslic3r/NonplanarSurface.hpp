#ifndef slic3r_NonplanarSurface_hpp_
#define slic3r_NonplanarSurface_hpp_

#include "libslic3r.h"
#include "NonplanarFacet.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ExPolygon.hpp"
#include "Geometry.hpp"
#include "ClipperUtils.hpp"

#include <cfenv>
#include <cmath>
#include <iostream>


namespace Slic3r {

typedef struct {
  float x;
  float y;
  float z;
} mesh_vertex;

typedef struct {
  mesh_vertex    max;
  mesh_vertex    min;
} mesh_stats;

class NonplanarSurface;

typedef std::vector<NonplanarSurface> NonplanarSurfaces;

class RayNP
{
public:
    Vec3f rayOrig;
    Vec3f rayEnd;      
};

class IntersectInfo 
{
public:
    IntersectInfo() 
    {
        bIntersect = false;
    }
    
    bool bIntersect;
    Vec3f intersectPos;
    float distance = -1.;
};


class RayTriangleIntersection
{
public:
    
    std::vector<IntersectInfo> checkMeshIntersection(std::vector<RayNP> p_rays, std::map<int, NonplanarFacet> p_faces)
    {
        faces.clear();
        rays.clear();
        intersecctInfos.clear();
        
        faces = p_faces;
        rays = p_rays;
        
        for (int j = 0; j < (int)rays.size(); j++)
        {
            RayNP ray = rays.at(j);
            Vec3f rayDir = ray.rayEnd;
            rayDir -= ray.rayOrig;
            
            for (int i = 0; i < (int)faces.size(); i++)
            {
                std::vector<Vec3f> verts;
                verts.push_back(faces.at(i).vertex[0]);
                verts.push_back(faces.at(i).vertex[1]);
                verts.push_back(faces.at(i).vertex[2]);        
                IntersectInfo intersecctInfo = judgeRayTri(ray.rayOrig, rayDir, verts);
                
                if (intersecctInfo.bIntersect)
                {
                    faces.at(i).bHit = true;
                    intersecctInfos.push_back(intersecctInfo);
                    break;
                }
            }
        }
        return intersecctInfos;
    }    

    std::vector<IntersectInfo>   intersecctInfos;
    std::vector<RayNP>             rays; 
    std::map<int, NonplanarFacet>  faces;
    
private:
    
    IntersectInfo judgeRayTri(Vec3f rayStart, Vec3f rayDir, std::vector<Vec3f> tri)
    {
        IntersectInfo result;   
        std::feclearexcept(FE_ALL_EXCEPT);
        rayDir.normalize();
        
        Vec3f triNorm = getNormal(tri);
        float vn = rayDir.dot(triNorm);
        
        Vec3f aa = rayStart - tri.at(0);
        float xpn = aa.dot(triNorm);
        float distance = -xpn / vn;
        
        if (distance < 0) return result; // behind ray origin. fail
        
        Vec3f hitPos = rayDir * distance + rayStart;
        
        Vec3f hit00 = hitPos - tri.at(0);
        Vec3f hit01 = tri.at(1) - tri.at(0);
        Vec3f cross0 = hit00.cross(hit01);
        if (cross0.dot(triNorm) > 0.000000001) return result;; // not in tri. fail
        
        Vec3f hit10 = hitPos - tri.at(1);
        Vec3f hit11 = tri.at(2) - tri.at(1);        
        Vec3f cross1 = hit10.cross(hit11);        
        if (cross1.dot(triNorm) > 0.000000001) return result;; // not in tri. fail        
        
        Vec3f hit20 = hitPos - tri.at(2);
        Vec3f hit21 = tri.at(0) - tri.at(2);        
        Vec3f cross2 = hit20.cross(hit21);        
        if (cross2.dot(triNorm) > 0.000000001) return result;; // not in tri. fail
        if (std::fetestexcept(FE_INVALID))
            return result;
        // intersect!
        result.bIntersect = true;
        result.intersectPos = hitPos;
        result.distance = distance;
        return result;
        
    }    
    
    Vec3f getNormal(std::vector<Vec3f> pverts)
    {
        Vec3f t0 = pverts.at(1)-pverts.at(0);
        Vec3f t1 = pverts.at(2)-pverts.at(0);
        Vec3f normal = t0.cross(t1);
        normal.normalize();
        return  normal;
    }
    
};

class NonplanarSurface
{
    public:
    std::map<int, NonplanarFacet> mesh;
    mesh_stats stats;
    NonplanarSurface() {};
    ~NonplanarSurface() {};
    NonplanarSurface(std::map<int, NonplanarFacet> &_mesh);
    bool operator==(const NonplanarSurface& other) const;
    void calculate_stats();
    void translate(float x, float y, float z);
    void scale(float factor);
    void scale(float versor[3]);
    void rotate_z(float angle);
    void debug_output();
    NonplanarSurfaces group_surfaces();
    void mark_neighbor_surfaces(int id);
    bool check_max_printing_height(float height);
    void check_printable_surfaces(float max_angle);
    bool check_surface_area();
    ExPolygons horizontal_projection() const;
    std::vector<IntersectInfo> intersect_rays(const Vec3f &src, const std::vector<Vec3f> &dirs);

};

};

#endif