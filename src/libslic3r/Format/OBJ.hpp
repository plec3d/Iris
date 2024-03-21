///|/ Copyright (c) Prusa Research 2017 - 2019 Tomáš Mészáros @tamasmeszaros, Vojtěch Bubník @bubnikv
///|/
///|/ ported from lib/Slic3r/Format/OBJ.pm:
///|/ Copyright (c) Prusa Research 2017 Vojtěch Bubník @bubnikv
///|/ Copyright (c) Slic3r 2012 - 2014 Alessandro Ranellucci @alranel
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_Format_OBJ_hpp_
#define slic3r_Format_OBJ_hpp_

#include <cstddef>
namespace Slic3r {

class TriangleMesh;
class Model;
class ModelObject;

// Load an OBJ file into a provided model.
extern bool load_obj(const char *path, std::vector<std::pair<std::string,TriangleMesh>> *mesh,  std::vector<std::pair<std::vector<int>,std::vector<Vec3f>>> *color_maps = nullptr);
extern bool load_obj(const char *path, Model *model, const char *object_name = nullptr);

extern bool store_obj(const char *path, TriangleMesh *mesh);
extern bool store_obj(const char *path, ModelObject *model);
extern bool store_obj(const char *path, Model *model);

}; // namespace Slic3r

#endif /* slic3r_Format_OBJ_hpp_ */
