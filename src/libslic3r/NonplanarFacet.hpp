#ifndef slic3r_NonplanarFacet_hpp_
#define slic3r_NonplanarFacet_hpp_

#include "libslic3r.h"
#include "Geometry.hpp"
#include "utility"

namespace Slic3r {

typedef struct {
  Vec3f    max;
  Vec3f    min;
} facet_stats;

class NonplanarFacet
{
    public:
    Vec3f vertex[3];
    Vec3f normal;
    int neighbor[3];
    facet_stats stats;
    bool marked = false;
    bool bHit;

    NonplanarFacet()
    {
        bHit = false;
    };
    ~NonplanarFacet() {};
    void calculate_stats();
    void translate(float x, float y, float z);
    void scale(float versor[3]);
    float calculate_surface_area();
    std::pair<float,Vec3f> find_closest_point(const Vec3f& x);
};

};

#endif