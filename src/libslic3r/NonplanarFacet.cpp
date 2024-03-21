#include "NonplanarFacet.hpp"
#include <utility>

namespace Slic3r {

void
NonplanarFacet::calculate_stats() {
    //calculate min and max values
    this->stats.min.x() = this->vertex[0].x();
    this->stats.min.y() = this->vertex[0].y();
    this->stats.min.z() = this->vertex[0].z();
    this->stats.max.x() = this->vertex[0].x();
    this->stats.max.y() = this->vertex[0].y();
    this->stats.max.z() = this->vertex[0].z();
    for(int j = 1; j < 3; j++) {
        this->stats.min.x() = std::min(this->stats.min.x(), this->vertex[j].x());
        this->stats.min.y() = std::min(this->stats.min.y(), this->vertex[j].y());
        this->stats.min.z() = std::min(this->stats.min.z(), this->vertex[j].z());
        this->stats.max.x() = std::max(this->stats.max.x(), this->vertex[j].x());
        this->stats.max.y() = std::max(this->stats.max.y(), this->vertex[j].y());
        this->stats.max.z() = std::max(this->stats.max.z(), this->vertex[j].z());
    }
}

void
NonplanarFacet::translate(float x, float y, float z)
{
    //translate facet
    for(int j = 0; j < 3; j++) {
      this->vertex[j].x() += x;
      this->vertex[j].y() += y;
      this->vertex[j].z() += z;
    }

    //translate min and max values
    this->stats.min.x() += x;
    this->stats.min.y() += y;
    this->stats.min.z() += z;
    this->stats.max.x() += x;
    this->stats.max.y() += y;
    this->stats.max.z() += z;
}

void
NonplanarFacet::scale(float versor[3])
{
    //scale facet
    for(int j = 0; j < 3; j++) {
        this->vertex[j].x() *= versor[0];
        this->vertex[j].y() *= versor[1];
        this->vertex[j].z() *= versor[2];
    }

    //scale min and max values
    this->stats.min.x() *= versor[0];
    this->stats.min.y() *= versor[1];
    this->stats.min.z() *= versor[2];
    this->stats.max.x() *= versor[0];
    this->stats.max.y() *= versor[1];
    this->stats.max.z() *= versor[2];
}

float
NonplanarFacet::calculate_surface_area()
{
    return Slic3r::Geometry::triangle_surface(
        Point(this->vertex[0].x(), this->vertex[0].y()), 
        Point(this->vertex[1].x(), this->vertex[1].y()), 
        Point(this->vertex[2].x(), this->vertex[2].y())
    );
}


std::pair<float,Vec3f> NonplanarFacet::find_closest_point(const Vec3f& x)
{
    Vec3f pt;
    Vec2f t;
	// source: real time collision detection
	// check if x in vertex region outside pa
	Vec3f ab = this->vertex[1] - this->vertex[0];
	Vec3f ac = this->vertex[2] - this->vertex[0];
	Vec3f ax = x - this->vertex[0];
	float d1 = ab.dot(ax);
	float d2 = ac.dot(ax);
	if (d1 <= 0.0f && d2 <= 0.0f) {
		// barycentric coordinates (1, 0, 0)
		t[0] = 1.0f;
		t[1] = 0.0f;
		pt = this->vertex[0];
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// check if x in vertex region outside pb
	Vec3f bx = x - this->vertex[1];
	float d3 = ab.dot(bx);
	float d4 = ac.dot(bx);
	if (d3 >= 0.0f && d4 <= d3) {
		// barycentric coordinates (0, 1, 0)
		t[0] = 0.0f;
		t[1] = 1.0f;
		pt = this->vertex[1];
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// check if x in vertex region outside pc
	Vec3f cx = x - this->vertex[2];
	float d5 = ab.dot(cx);
	float d6 = ac.dot(cx);
	if (d6 >= 0.0f && d5 <= d6) {
		// barycentric coordinates (0, 0, 1)
		t[0] = 0.0f;
		t[1] = 0.0f;
		pt = this->vertex[2];
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// check if x in edge region of ab, if so return projection of x onto ab
	float vc = d1*d4 - d3*d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
		// barycentric coordinates (1 - v, v, 0)
		float v = d1/(d1 - d3);
		t[0] = 1.0f - v;
		t[1] = v;
		pt = this->vertex[0] + ab*v;
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// check if x in edge region of ac, if so return projection of x onto ac
	float vb = d5*d2 - d1*d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		// barycentric coordinates (1 - w, 0, w)
		float w = d2/(d2 - d6);
		t[0] = 1.0f - w;
		t[1] = 0.0f;
		pt = this->vertex[0] + ac*w;
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// check if x in edge region of bc, if so return projection of x onto bc
	float va = d3*d6 - d5*d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
		// barycentric coordinates (0, 1 - w, w)
		float w = (d4 - d3)/((d4 - d3) + (d5 - d6));
		t[0] = 0.0f;
		t[1] = 1.0f - w;
		pt = this->vertex[1] + (this->vertex[2]- this->vertex[1])*w;
		return std::pair<float,Vec3f>{(x - pt).norm(),pt};
	}

	// x inside face region. Compute pt through its barycentric coordinates (u, v, w)
	float denom = 1.0f/(va + vb + vc);
	float v = vb*denom;
	float w = vc*denom;
	t[0] = 1.0f - v - w;
	t[1] = v;

	pt = this->vertex[0] + ab*v + ac*w; //= u*a + v*b + w*c, u = va*denom = 1.0f - v - w
	return std::pair<float,Vec3f>{(x - pt).norm(),pt};
}

}
