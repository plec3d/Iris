#include "FillSpiral.hpp"
#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Polygon.hpp"
#include "../Polyline.hpp"
#include "../Surface.hpp"

namespace Slic3r {

void FillSpiral::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction,
    const Polyline                  &pedestal,
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    BoundingBox bbox = expolygon.contour.bounding_box();
    coord_t min_spacing = scale_(this->spacing);
    coord_t distance = coord_t(min_spacing / params.density);
    double bridge_angl = (int(direction.first*100)%int(2*M_PI))/100;
    double center_x = bbox.min(0) + bbox.size()(0)/2;//*std::cos(direction.first);//* bridge_angl;// * angleSignTable[bridge_sign][0];//*std::abs((angleSignTable[bridge_sign][0]>0?1:0)-bridge_angl);*/
    double center_y = bbox.min(1) + bbox.size()(1)/2;//*std::sin(direction.first);//* bridge_angl;//+ radius;// * angleSignTable[bridge_sign][1];//*std::abs((angleSignTable[bridge_sign][0]>0?0:1)-bridge_angl);*/
    double angle = direction.first-.5*M_PI;
        center_x = pedestal.bounding_box().center().x();
        center_y = pedestal.bounding_box().center().y();
        // Now unwind the spiral.
        Pointfs out;
        Polylines polylines;
        //FIXME Vojtech: If used as a solid infill, there is a gap left at the center.
        // Radius to achieve.
        coord_t min_x = coord_t(ceil(coordf_t(bbox.min.x()) / distance));
        coord_t min_y = coord_t(ceil(coordf_t(bbox.min.y()) / distance));
        coord_t max_x = coord_t(ceil(coordf_t(bbox.max.x()) / distance));
        coord_t max_y = coord_t(ceil(coordf_t(bbox.max.y()) / distance));
        coordf_t rmax = std::sqrt(coordf_t(max_x)*coordf_t(max_x)+coordf_t(max_y)*coordf_t(max_y)) * std::sqrt(3.) + 1.5;
        // Now unwind the spiral.
        //coordf_t a = 1.;
        coordf_t b = 1./(2.*M_PI);
        coordf_t theta = 0.;
        coordf_t r = 1;
	    out.emplace_back(0,0);
        //FIXME Vojtech: If used as a solid infill, there is a gap left at the center.
        while (r < rmax) {
            // Discretization angle to achieve a discretization error lower than resolution.
            theta += 2. * acos(1. - params.resolution / r);
            r = /*a + */b * theta;
            out.emplace_back(r * cos(theta), r * sin(theta));
        }
        if (out.size() >= 2) {
            // Convert points to a polyline, upscale.
            Polylines polylines(1, Polyline());
            Polyline &polyline = polylines.front();
            polyline.points.reserve(out.size());
            for (const Vec2d &pt : out)
                polyline.points.emplace_back(
                    coord_t(floor(center_x + pt.x() * distance + 0.5)),
                    coord_t(floor(center_y + pt.y() * distance + 0.5)));
            polylines = intersection_pl(polylines, expolygon);
	    std::sort(polylines.begin(), polylines.end(),
              [center_x,center_y](const auto &i1, const auto &i2) {
                  return i1.distance_to(Point(center_x,center_y)) < i2.distance_to(Point(center_x,center_y));/* .at(0).x() >
                             i2.at(0).x() &&
                         i1.at(0).y() > i2.at(0).y();*/
              });
            append(polylines_out, std::move(polylines));
        }
    }
} // namespace Slic3r