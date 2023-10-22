#include "FillArc.hpp"
#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Polygon.hpp"
#include "../Polyline.hpp"
#include "../Surface.hpp"

namespace Slic3r {

void FillArc::_fill_surface_single(
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
    coord_t radius = std::max(bbox.max(0)-bbox.min(0),bbox.max(1)-bbox.min(1))*params.config->arc_radius;
    double center_x = bbox.min(0) + bbox.size()(0)/2;//*std::cos(direction.first);//* bridge_angl;// * angleSignTable[bridge_sign][0];//*std::abs((angleSignTable[bridge_sign][0]>0?1:0)-bridge_angl);*/
    double center_y = bbox.min(1) + bbox.size()(1)/2;//*std::sin(direction.first);//* bridge_angl;//+ radius;// * angleSignTable[bridge_sign][1];//*std::abs((angleSignTable[bridge_sign][0]>0?0:1)-bridge_angl);*/
    // do a raytrace to find the place where the angle from center hits the contour and translate the arc's center to that place
    Polyline intersect_line;
    Polylines intersect_lines;
    ClipperLib::Clipper clipper;
    Vec2d pt;
    double angle = direction.first-.5*M_PI;
    intersect_line.points.emplace_back((center_x),(center_y));
    intersect_line.points.emplace_back(((double)(center_x + bbox.size().x()*double(params.config->arc_infill_raylen.value)*std::cos(1*angle)))
                            ,((double)(center_y + bbox.size().y()*double(params.config->arc_infill_raylen.value)*std::sin(1*angle))));
    intersect_line.points.emplace_back(((double)(center_x + bbox.size().x()*(double(params.config->arc_infill_raylen.value)+.1)*std::cos(1*angle)))
                            ,((double)(center_y + bbox.size().y()*(double(params.config->arc_infill_raylen.value)+.1)*std::sin(1*angle))));
    intersect_line.points.emplace_back((center_x),(center_y));
    intersect_lines = intersection_pl(intersect_line,expolygon);
    if(!intersect_lines.empty()){
        center_x = (direction.first<M_PI||direction.first>2*M_PI)?(std::max(std::max(intersect_lines.front().back().x(),intersect_lines.front().front().x()),std::max(intersect_lines.back().back().x(),intersect_lines.back().front().x()))):
                    (std::min(std::min(intersect_lines.front().back().x(),intersect_lines.front().front().x()),std::min(intersect_lines.back().back().x(),intersect_lines.back().front().x())));
        center_y = (direction.first<1.5*M_PI)?(std::max(std::max(intersect_lines.front().back().y(),intersect_lines.front().front().y()),std::max(intersect_lines.back().back().y(),intersect_lines.back().front().y()))):
                    (std::min(std::min(intersect_lines.front().back().y(),intersect_lines.front().front().y()),std::min(intersect_lines.back().back().y(),intersect_lines.back().front().y())));
    }else{
    intersect_line.clear();
    intersect_line.points.emplace_back((center_x),(center_y));
    intersect_line.points.emplace_back(((double)(bbox.size().x()*double(params.config->arc_infill_raylen.value)*std::cos(-1*angle)))
                            ,((double)(bbox.size().y()*double(params.config->arc_infill_raylen.value)*std::sin(-1*angle))));
    intersect_line.points.emplace_back(((double)(bbox.size().x()*(double(params.config->arc_infill_raylen.value)+.1)*std::cos(-1*angle)))
                            ,((double)(bbox.size().y()*(double(params.config->arc_infill_raylen.value)+.1)*std::sin(-1*angle))));
    intersect_line.points.emplace_back((center_x),(center_y));
    intersect_lines = intersection_pl(intersect_line,expolygon);
    if(!intersect_lines.empty()){
        center_x = (direction.first<M_PI||direction.first>2*M_PI)?(std::max(std::max(intersect_lines.front().back().x(),intersect_lines.front().front().x()),std::max(intersect_lines.back().back().x(),intersect_lines.back().front().x()))):
                    (std::min(std::min(intersect_lines.front().back().x(),intersect_lines.front().front().x()),std::min(intersect_lines.back().back().x(),intersect_lines.back().front().x())));
        center_y = (direction.first<1.5*M_PI)?(std::max(std::max(intersect_lines.front().back().y(),intersect_lines.front().front().y()),std::max(intersect_lines.back().back().y(),intersect_lines.back().front().y()))):
                    (std::min(std::min(intersect_lines.front().back().y(),intersect_lines.front().front().y()),std::min(intersect_lines.back().back().y(),intersect_lines.back().front().y())));
    }else{
        center_x = bbox.center().x()+bbox.size().x()/2*std::cos((angle+(bridge_angl)*M_PI));
        center_y = bbox.center().y()+bbox.size().y()/2*std::sin((angle+(bridge_angl)*M_PI));
        }
    }
    ExPolygon expolygons;
    for (double angle=0.; angle<=2.05*M_PI; angle+=0.05*M_PI)
        expolygons.contour.points.push_back(Point((double)center_x + radius * cos(angle),
                                            (double)center_y + radius * sin(angle)));
    Polygons  loops = to_polygons(expolygons);
    ExPolygons  last { std::move(expolygons) };    
    Polyline line;
    Polyline line0;
    Polylines lines;
    Polylines polylines;
    // make the new outlines of the circles
    while (!last.empty()) {
        int i = 0;
        last = offset2_ex(last, -(distance + min_spacing/2), +min_spacing/2);//,ClipperLib::jtRound);  
        
        for (ExPolygons::const_iterator it = last.begin(); it != last.end();
            ++it) {
        line.clear();
        line0.clear();
        for (const Point &point : it->contour.points) {
            pt = point.cast<double>();
            // add_offset_line(pt);
            /*pt += Vec2d(0.5 - (pt.x() < 0),
                        0.5 - (pt.y() < 0));*/
            if(!(!expolygon.holes.empty() &&(expolygon.contains_h(direction.second)))
                && i++<it->contour.points.size()*std::tan((direction.first-.25*M_PI)/6)) line0.points.push_back(Point(pt.x(),pt.y()));
            else
        line.points.push_back(Point(pt.x(),pt.y()));/*ClipperLib::cInt(pt.x()),
                                ClipperLib::cInt(pt.y()));
        */}
        if(!line0.empty()) for (const Point &point : line0) {
            pt = point.cast<double>();
            // add_offset_line(pt);
            //pt += Vec2d(0.5 - (pt.x() < 0),
            //            0.5 - (pt.y() < 0));
            line.points.push_back(Point(pt.x(),pt.y()));
        }
        if(line.points.front()!=line.points.back()) line.points.push_back(line.points.front());
        lines.push_back(line);
        }
    }
    polylines = intersection_pl(lines, expolygon);
    std::sort(polylines.begin(), polylines.end(),
            [center_x,center_y](const auto &i1, const auto &i2) {
        return i1.distance_from_edges(Point(center_x,center_y)) < i2.distance_from_edges(Point(center_x,center_y));
            });
    Polylines chained;
    Point ptc(center_x,center_y);
    chained = chain_polylines(std::move(polylines),&ptc);
    append(polylines_out, std::move(polylines.size()==chained.size()?chained:polylines));
}
} // namespace Slic3r