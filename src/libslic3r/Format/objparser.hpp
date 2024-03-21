///|/ Copyright (c) Prusa Research 2017 - 2019 Tomáš Mészáros @tamasmeszaros, Vojtěch Bubník @bubnikv
///|/ Copyright (c) 2017 Joseph Lenox @lordofhyphens
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_Format_objparser_hpp_
#define slic3r_Format_objparser_hpp_

#include <string>
#include <vector>
#include <istream>

namespace ObjParser {

struct ObjUseMtl
{
	int			vertexIdxFirst;
	int			normalIdxFirst;
	int			textureIdxFirst;

	std::string name;
};

inline bool operator==(const ObjUseMtl &v1, const ObjUseMtl &v2)
{
	return 
		v1.vertexIdxFirst	== v2.vertexIdxFirst	&&
		v1.normalIdxFirst	== v2.normalIdxFirst	&&
		v1.textureIdxFirst	== v2.textureIdxFirst	&&
		v1.name.compare(v2.name) == 0;
}

struct ObjObject
{
	int			coordIdxFirst;
	int			vertexIdxFirst;
	int			normalIdxFirst;
	int			textureIdxFirst;
	std::string name;
};

inline bool operator==(const ObjObject &v1, const ObjObject &v2)
{
	return 
		v1.vertexIdxFirst	== v2.vertexIdxFirst	&&
		v1.normalIdxFirst	== v2.normalIdxFirst	&&
		v1.textureIdxFirst	== v2.textureIdxFirst	&& 
		v1.name.compare(v2.name) == 0;
}

struct ObjGroup
{
	int			coordIdxFirst;
	int			vertexIdxFirst;
	int			normalIdxFirst;
	int			textureIdxFirst;
	std::string name;
};

inline bool operator==(const ObjGroup &v1, const ObjGroup &v2)
{
	return 
		v1.coordIdxFirst	== v2.coordIdxFirst	&& 
		v1.vertexIdxFirst	== v2.vertexIdxFirst	&&
		v1.normalIdxFirst	== v2.normalIdxFirst	&&
		v1.textureIdxFirst	== v2.textureIdxFirst	&&
		v1.name.compare(v2.name) == 0;
}

struct ObjSmoothingGroup
{
	int			vertexIdxFirst;
	int			normalIdxFirst;
	int			textureIdxFirst;
	int			smoothingGroupID;
};

inline bool operator==(const ObjSmoothingGroup &v1, const ObjSmoothingGroup &v2)
{
	return 
		v1.vertexIdxFirst	== v2.vertexIdxFirst	&&
		v1.normalIdxFirst	== v2.normalIdxFirst	&&
		v1.textureIdxFirst	== v2.textureIdxFirst	&&
		v1.smoothingGroupID == v2.smoothingGroupID;
}

struct ObjData {
	// Version of the data structure for load / store in the private binary format.
	int								version;
	// x, y, z
	std::vector<Slic3r::Vec3f>				coordinates;
	// x, y, z
	std::vector<Slic3r::Vec3f>				colors;
	// u, v, w
	std::vector<Slic3r::Vec3f>				textureCoordinates;
	// x, y, z
	std::vector<Slic3r::Vec3f>				normals;
	// u, v, w
	std::vector<Slic3r::Vec3f>				parameters;

	std::vector<std::string>		mtllibs;
	std::vector<ObjUseMtl>			usemtls;
	std::vector<ObjObject>			objects;
	std::vector<ObjGroup>			groups;
	std::vector<ObjSmoothingGroup>	smoothingGroups;
	std::vector<bool>				is_quad;
	// List of faces, delimited by an ObjVertex with all members set to -1.
	std::vector<Slic3r::Vec3i>			faces;
	std::vector<Slic3r::Vec3i>			face_normals;
	std::vector<Slic3r::Vec3i>			face_texts;
};

extern bool objparse(const char *path, ObjData &data);
extern bool objparse(std::istream &stream, ObjData &data);

extern bool objbinsave(const char *path, const ObjData &data);

extern bool objbinload(const char *path, ObjData &data);

extern bool objequal(const ObjData &data1, const ObjData &data2);

} // namespace ObjParser

#endif /* slic3r_Format_objparser_hpp_ */
