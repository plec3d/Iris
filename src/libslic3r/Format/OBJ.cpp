///|/ Copyright (c) Prusa Research 2017 - 2021 Enrico Turri @enricoturri1966, Vojtěch Bubník @bubnikv, Tomáš Mészáros @tamasmeszaros
///|/
///|/ ported from lib/Slic3r/Format/OBJ.pm:
///|/ Copyright (c) Prusa Research 2017 Vojtěch Bubník @bubnikv
///|/ Copyright (c) Slic3r 2012 - 2014 Alessandro Ranellucci @alranel
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "../libslic3r.h"
#include "../Model.hpp"
#include "../Color.hpp"
#include "../TriangleMesh.hpp"

#include "OBJ.hpp"
#include "objparser.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <boost/log/trivial.hpp>

#ifdef _WIN32
#define DIR_SEPARATOR '\\'
#else
#define DIR_SEPARATOR '/'
#endif

namespace Slic3r {

bool load_obj(const char *path, std::vector<std::pair<std::string,TriangleMesh>> *meshptr, std::vector<std::pair<std::vector<int>,std::vector<Vec3f>>> *color_maps)
{
    if (meshptr == nullptr)
        return false;
    
    // Parse the OBJ file.
    ObjParser::ObjData data;
    if (! ObjParser::objparse(path, data)) {
        BOOST_LOG_TRIVIAL(error) << "load_obj: failed to parse " << path;
        return false;
    }
    int i = 0;
    for(auto &group: data.groups){

        std::pair<std::vector<int>,std::vector<Vec3f>> color_map;
        int nextVertexId = i<data.groups.size()-1?data.groups[++i].vertexIdxFirst:data.faces.size();
        int nextCoordId = i<data.groups.size()-1?data.groups[i].coordIdxFirst:data.coordinates.size();
        //if(data.is_quad.back())
        //    nextVertexId--;
        // Count the faces and verify, that all faces are triangular.
        size_t num_faces = nextVertexId - group.vertexIdxFirst;
        size_t num_quads = 0;
        
        // Convert ObjData into indexed triangle set.
        indexed_triangle_set its;
        its.vertices.reserve(nextCoordId - group.coordIdxFirst);
        its.indices.reserve(num_faces + num_quads);
        for (size_t j = group.coordIdxFirst; j < nextCoordId; j++) {
            its.vertices.emplace_back(std::move(data.coordinates[j]));
            
            if(color_maps != nullptr && data.colors.size() > j && data.colors[j][0] >= 0.f && data.colors[j][1] >= 0.f && data.colors[j][2] >= 0.f){
                std::vector<Vec3f>::iterator it = std::find_if (color_map.second.begin(), color_map.second.end(), [&data, j](Vec3f i) {
                    return (i == data.colors[j] || (i - data.colors[j]).norm() < 0.3);
                });
                if(it!=color_map.second.end())
                    color_map.first.push_back(std::distance(color_map.second.begin(), it));
                else {                
                    color_map.second.push_back(std::move(data.colors[j]));
                    color_map.first.push_back(color_map.second.size()-1);
                    BOOST_LOG_TRIVIAL(info) << "Colors mapped: " << color_map.second.back()[0] << " "<< color_map.second.back()[1]<< " "<< color_map.second.back()[2];
                }
            }
        }
         BOOST_LOG_TRIVIAL(info) << "Colors mapped: " << color_map.second.size();
        // fill the indices
        for (int i = group.vertexIdxFirst; i < nextVertexId; i++){
            //BOOST_LOG_TRIVIAL(info) << data.faces.size()-1 << " " << data.faces[i][0];
            if(data.faces[i][0] >= 0 && data.faces[i][1] >= 0 && data.faces[i][2] >= 0)
                its.indices.emplace_back(std::move(data.faces[i]));
            //if(data.is_quad.at(i))
            //    its.indices.emplace_back(std::move(data.faces[++i]));
        }
        //BOOST_LOG_TRIVIAL(info) << "Colors mapped: " << color_map.second.size();
        if (!its.indices.empty())
            meshptr->push_back(std::make_pair(group.name,TriangleMesh(std::move(its))));
        else
            continue;
        if (meshptr->back().second.empty()) {
            //BOOST_LOG_TRIVIAL(error) << "load_obj: This OBJ file couldn't be read because it's empty. " << path;
            return false;
        }
        //BOOST_LOG_TRIVIAL(info) << "Colors mapped: " << color_map.second.size();
        if (meshptr->back().second.volume() < 0)
            meshptr->back().second.flip_triangles();
        color_maps->push_back(color_map);
    }
    return true;
}

bool load_obj(const char *path, Model *model, const char *object_name_in)
{
    std::vector<std::pair<std::string,TriangleMesh>> meshes;
    std::vector<std::pair<std::vector<int>,std::vector<Vec3f>>> color_maps;

    bool ret = load_obj(path, &meshes, &color_maps);
    
    if (ret) {
        int j = 0;
        for(auto &mesh: meshes){
            std::pair<std::vector<int>,std::vector<Vec3f>> &color_map = color_maps.at(j++);
            std::string  object_name = mesh.first;
            /*if (object_name_in == nullptr) {
                const char *last_slash = strrchr(path, DIR_SEPARATOR);
                object_name.assign((last_slash == nullptr) ? path : last_slash + 1);
            } else
            object_name.assign(object_name_in);*/
        
            model->add_object(object_name.c_str(), path, std::move(mesh.second));
            if(color_map.first.size()>0){
                int i = 0;
                for(Vec3f &col: color_map.second){
                    std::vector<std::pair<int,std::string>>::iterator it = std::find_if (model->objects.back()->mmu_colors.begin(), model->objects.back()->mmu_colors.end(), [&col](std::pair<int,std::string> i) {
                        ColorRGB color;
                        decode_color(i.second,color);
                        Vec3f colo(color.r(),color.g(),color.b());
                        return (colo == col || (colo - col).norm() < .5);
                    });
                    if(it==model->objects.back()->mmu_colors.end())
                        model->objects.back()->mmu_colors.push_back(std::make_pair(-1,encode_color(ColorRGB((float)col[0],(float)col[1],(float)col[2]))));
                    
                    for(int col: color_map.first){
                        if(col == i)
                            col =  it==model->objects.back()->mmu_colors.end()? model->objects.back()->mmu_colors.size()-1:
                                std::distance(model->objects.back()->mmu_colors.begin(), it);
                    }
                    i++;
                }
                
                model->objects.back()->mmu_indices_map = std::move(color_map.first);
                BOOST_LOG_TRIVIAL(info) << "Colors mapped: " << model->objects.back()->mmu_indices_map.size();
            }
        }
    }
    
    return ret;
}

bool store_obj(const char *path, TriangleMesh *mesh)
{
    //FIXME returning false even if write failed.
    mesh->WriteOBJFile(path);
    return true;
}

bool store_obj(const char *path, ModelObject *model_object)
{
    TriangleMesh mesh = model_object->mesh();
    return store_obj(path, &mesh);
}

bool store_obj(const char *path, Model *model)
{
    TriangleMesh mesh = model->mesh();
    return store_obj(path, &mesh);
}

}; // namespace Slic3r
