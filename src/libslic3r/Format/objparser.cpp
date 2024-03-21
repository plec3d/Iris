///|/ Copyright (c) Prusa Research 2017 - 2021 Lukáš Matěna @lukasmatena, Vojtěch Bubník @bubnikv, Tomáš Mészáros @tamasmeszaros, Enrico Turri @enricoturri1966
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/

#include "objparser.hpp"

namespace ObjParser {

std::vector<std::string> split_string(std::string s, std::string delimiter) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
		token = s.substr (pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		if(token != "") res.push_back (token);
	}

	res.push_back (s.substr (pos_start));
	return res;
}

bool objparse(const char *path, ObjData &data)
{
	std::ifstream fin(path);
	std::string str;
	std::vector<std::string> words;
	while (!fin.eof()) {
		std::getline(fin, str);
		words = split_string(str, " ");

		if (words.size() == 4 && words[0] == "v") {
			Slic3r::Vec3f v_tmp(
				atof(words[1].c_str()),
				atof(words[2].c_str()),
				atof(words[3].c_str()));
			if(std::isnan(v_tmp[0]))
				v_tmp[0] = data.coordinates.back()[0];
			if(std::isnan(v_tmp[1]))
				v_tmp[1] = data.coordinates.back()[1];
			if(std::isnan(v_tmp[2]))
				v_tmp[2] = data.coordinates.back()[2];
			data.coordinates.push_back(std::move(v_tmp));
		}else if (words.size() == 7 && words[0] == "v") {
			Slic3r::Vec3f v_tmp(
				atof(words[1].c_str()),
				atof(words[2].c_str()),
				atof(words[3].c_str()));
			if(std::isnan(v_tmp[0]))
				v_tmp[0] = data.coordinates.back()[0];
			if(std::isnan(v_tmp[1]))
				v_tmp[1] = data.coordinates.back()[1];
			if(std::isnan(v_tmp[2]))
				v_tmp[2] = data.coordinates.back()[2];
			data.coordinates.push_back(std::move(v_tmp));
			data.colors.push_back(Slic3r::Vec3f(
				atof(words[4].c_str()),
				atof(words[5].c_str()),
				atof(words[6].c_str()))
			);
		}
		/*else if (words.size() == 3 && words[0] == "vt") {
			data.textureCoordinates.push_back(Slic3r::Vec3f(
				atof(words[1].c_str()),
				atof(words[2].c_str()))
			);
		}*/
		else if (words.size() == 4 && words[0] == "vn") {
			data.normals.push_back(Slic3r::Vec3f(
				atof(words[1].c_str()),
				atof(words[2].c_str()),
				atof(words[3].c_str()))
			);
		}
		else if (words.size() == 2 && words[0] == "g") {
			ObjGroup group;
			group.coordIdxFirst = (int)data.coordinates.size();
			group.vertexIdxFirst = (int)data.faces.size();
			group.normalIdxFirst = (int)data.normals.size();
			group.textureIdxFirst = (int)data.textureCoordinates.size();
			group.name = words[1];
			data.groups.push_back(group);
		}
		else if (words.size() == 2 && words[0] == "o") {
			ObjObject object;
			object.coordIdxFirst = (int)data.coordinates.size();
			object.vertexIdxFirst = (int)data.faces.size();
			object.normalIdxFirst = (int)data.normals.size();
			object.textureIdxFirst = (int)data.textureCoordinates.size();
			object.name = words[1];
			data.objects.push_back(object);
		}
		else if (words.size() >= 4 && words[0] == "f") {
			Slic3r::Vec3i face_tmp(-1,-1,-1);
			Slic3r::Vec3i normal_tmp(-1,-1,-1);
			Slic3r::Vec3i text_tmp(-1,-1,-1);
			for (int i = 1; i < words.size(); ++i) {
				std::vector str_tmp = split_string(words[i], "/");
				if(i==4){
					face_tmp[1] = face_tmp[2];
					normal_tmp[1] = normal_tmp[2];
					text_tmp[1] = text_tmp[2];
					face_tmp[2] = atoi(str_tmp[0].c_str()) - 1 - data.groups.back().coordIdxFirst;
					normal_tmp[2] = atoi(str_tmp[1].c_str()) - 1 - data.groups.back().normalIdxFirst;
					text_tmp[2] = atoi(str_tmp[2].c_str()) - 1 - data.groups.back().textureIdxFirst;
					data.is_quad.at(data.faces.size()-1) = true;
					data.is_quad.push_back(true);
					data.faces.push_back(std::move(face_tmp));
					data.face_normals.push_back(std::move(normal_tmp));
					data.face_texts.push_back(std::move(text_tmp));
				}else{
					face_tmp[i-1] = atoi(str_tmp[0].c_str()) - 1 - data.groups.back().coordIdxFirst;
					normal_tmp[i-1] = atoi(str_tmp[1].c_str()) - 1 - data.groups.back().normalIdxFirst;
					text_tmp[i-1] = atoi(str_tmp[2].c_str()) - 1 - data.groups.back().textureIdxFirst;
					data.faces.push_back(std::move(face_tmp));
					data.face_normals.push_back(std::move(normal_tmp));
					data.face_texts.push_back(std::move(text_tmp));
					data.is_quad.push_back(false);
				}
			}
		}
	}
	return true;
}

} // namespace ObjParser
