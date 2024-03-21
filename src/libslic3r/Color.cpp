///|/ Copyright (c) Prusa Research 2021 - 2022 Vojtěch Bubník @bubnikv, Enrico Turri @enricoturri1966
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "libslic3r.h"
#include "Color.hpp"
#include <boost/log/trivial.hpp>

#include <random>
#include <algorithm>
#include <vector>

static const float INV_255 = 1.0f / 255.0f;

namespace Slic3r {

// Conversion from RGB to HSV color space
// The input RGB values are in the range [0, 1]
// The output HSV values are in the ranges h = [0, 360], and s, v = [0, 1]
static void RGBtoHSV(float r, float g, float b, float& h, float& s, float& v)
{
	assert(0.0f <= r && r <= 1.0f);
	assert(0.0f <= g && g <= 1.0f);
	assert(0.0f <= b && b <= 1.0f);

	const float max_comp = std::max(std::max(r, g), b);
	const float min_comp = std::min(std::min(r, g), b);
	const float delta = max_comp - min_comp;

	if (delta > 0.0f) {
		if (max_comp == r)
			h = 60.0f * (std::fmod(((g - b) / delta), 6.0f));
		else if (max_comp == g)
			h = 60.0f * (((b - r) / delta) + 2.0f);
		else  // max_comp == b
			h = 60.0f * (((r - g) / delta) + 4.0f);

		s = (max_comp > 0.0f) ? delta / max_comp : 0.0f;
	}
	else {
		h = 0.0f;
		s = 0.0f;
	}
	v = max_comp;

	while (h < 0.0f) { h += 360.0f; }
	while (h > 360.0f) { h -= 360.0f; }

	assert(0.0f <= s && s <= 1.0f);
	assert(0.0f <= v && v <= 1.0f);
	assert(0.0f <= h && h <= 360.0f);
}

// Conversion from HSV to RGB color space
// The input HSV values are in the ranges h = [0, 360], and s, v = [0, 1]
// The output RGB values are in the range [0, 1]
static void HSVtoRGB(float h, float s, float v, float& r, float& g, float& b)
{
	assert(0.0f <= s && s <= 1.0f);
	assert(0.0f <= v && v <= 1.0f);
	assert(0.0f <= h && h <= 360.0f);

	const float chroma = v * s;
	const float h_prime = std::fmod(h / 60.0f, 6.0f);
	const float x = chroma * (1.0f - std::abs(std::fmod(h_prime, 2.0f) - 1.0f));
	const float m = v - chroma;

	if (0.0f <= h_prime && h_prime < 1.0f) {
		r = chroma;
		g = x;
		b = 0.0f;
	}
	else if (1.0f <= h_prime && h_prime < 2.0f) {
		r = x;
		g = chroma;
		b = 0.0f;
	}
	else if (2.0f <= h_prime && h_prime < 3.0f) {
		r = 0.0f;
		g = chroma;
		b = x;
	}
	else if (3.0f <= h_prime && h_prime < 4.0f) {
		r = 0.0f;
		g = x;
		b = chroma;
	}
	else if (4.0f <= h_prime && h_prime < 5.0f) {
		r = x;
		g = 0.0f;
		b = chroma;
	}
	else if (5.0f <= h_prime && h_prime < 6.0f) {
		r = chroma;
		g = 0.0f;
		b = x;
	}
	else {
		r = 0.0f;
		g = 0.0f;
		b = 0.0f;
	}

	r += m;
	g += m;
	b += m;

	assert(0.0f <= r && r <= 1.0f);
	assert(0.0f <= g && g <= 1.0f);
	assert(0.0f <= b && b <= 1.0f);
}

class Randomizer
{
	std::random_device m_rd;

public:
	float random_float(float min, float max) {
		std::mt19937 rand_generator(m_rd());
		std::uniform_real_distribution<float> distrib(min, max);
		return distrib(rand_generator);
	}
};

ColorRGB::ColorRGB(float r, float g, float b)
: m_data({ std::clamp(r, 0.0f, 1.0f), std::clamp(g, 0.0f, 1.0f), std::clamp(b, 0.0f, 1.0f) })
{
}

ColorRGB::ColorRGB(unsigned char r, unsigned char g, unsigned char b)
: m_data({ std::clamp(r * INV_255, 0.0f, 1.0f), std::clamp(g * INV_255, 0.0f, 1.0f), std::clamp(b * INV_255, 0.0f, 1.0f) })
{
}

bool ColorRGB::operator < (const ColorRGB& other) const
{
	for (size_t i = 0; i < 3; ++i) {
		if (m_data[i] < other.m_data[i])
			return true;
		else if (m_data[i] > other.m_data[i])
			return false;
	}

	return false;
}

bool ColorRGB::operator > (const ColorRGB& other) const
{
	for (size_t i = 0; i < 3; ++i) {
		if (m_data[i] > other.m_data[i])
			return true;
		else if (m_data[i] < other.m_data[i])
			return false;
	}

	return false;
}

ColorRGB ColorRGB::operator + (const ColorRGB& other) const
{
	ColorRGB ret;
	for (size_t i = 0; i < 3; ++i) {
		ret.m_data[i] = std::clamp(m_data[i] + other.m_data[i], 0.0f, 1.0f);
	}
	return ret;
}

ColorRGB ColorRGB::operator * (float value) const
{
	assert(value >= 0.0f);
	ColorRGB ret;
	for (size_t i = 0; i < 3; ++i) {
		ret.m_data[i] = std::clamp(value * m_data[i], 0.0f, 1.0f);
	}
	return ret;
}

ColorRGBA::ColorRGBA(float r, float g, float b, float a)
: m_data({ std::clamp(r, 0.0f, 1.0f), std::clamp(g, 0.0f, 1.0f), std::clamp(b, 0.0f, 1.0f), std::clamp(a, 0.0f, 1.0f) })
{
}

ColorRGBA::ColorRGBA(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
: m_data({ std::clamp(r * INV_255, 0.0f, 1.0f), std::clamp(g * INV_255, 0.0f, 1.0f), std::clamp(b * INV_255, 0.0f, 1.0f), std::clamp(a * INV_255, 0.0f, 1.0f) })
{
}

bool ColorRGBA::operator < (const ColorRGBA& other) const
{
	for (size_t i = 0; i < 3; ++i) {
		if (m_data[i] < other.m_data[i])
			return true;
		else if (m_data[i] > other.m_data[i])
			return false;
	}

	return false;
}

bool ColorRGBA::operator > (const ColorRGBA& other) const
{
	for (size_t i = 0; i < 3; ++i) {
		if (m_data[i] > other.m_data[i])
			return true;
		else if (m_data[i] < other.m_data[i])
			return false;
	}

	return false;
}

ColorRGBA ColorRGBA::operator + (const ColorRGBA& other) const
{
	ColorRGBA ret;
	for (size_t i = 0; i < 3; ++i) {
		ret.m_data[i] = std::clamp(m_data[i] + other.m_data[i], 0.0f, 1.0f);
	}
	return ret;
}

ColorRGBA ColorRGBA::operator * (float value) const
{
	assert(value >= 0.0f);
	ColorRGBA ret;
	for (size_t i = 0; i < 3; ++i) {
		ret.m_data[i] = std::clamp(value * m_data[i], 0.0f, 1.0f);
	}
	ret.m_data[3] = m_data[3];
	return ret;
}

ColorRGB operator * (float value, const ColorRGB& other) { return other * value; }
ColorRGBA operator * (float value, const ColorRGBA& other) { return other * value; }

ColorRGB lerp(const ColorRGB& a, const ColorRGB& b, float t)
{
	assert(0.0f <= t && t <= 1.0f);
	return (1.0f - t) * a + t * b;
}

ColorRGBA lerp(const ColorRGBA& a, const ColorRGBA& b, float t)
{
	assert(0.0f <= t && t <= 1.0f);
	return (1.0f - t) * a + t * b;
}

ColorRGB complementary(const ColorRGB& color)
{
	return { 1.0f - color.r(), 1.0f - color.g(), 1.0f - color.b() };
}

ColorRGBA complementary(const ColorRGBA& color)
{
	return { 1.0f - color.r(), 1.0f - color.g(), 1.0f - color.b(), color.a() };
}

ColorRGB saturate(const ColorRGB& color, float factor)
{
	float h, s, v;
	RGBtoHSV(color.r(), color.g(), color.b(), h, s, v);
	s = std::clamp(s * factor, 0.0f, 1.0f);
	float r, g, b;
	HSVtoRGB(h, s, v, r, g, b);
	return { r, g, b };
}

ColorRGBA saturate(const ColorRGBA& color, float factor)
{
	return to_rgba(saturate(to_rgb(color), factor), color.a());
}

ColorRGB opposite(const ColorRGB& color)
{
	float h, s, v;
	RGBtoHSV(color.r(), color.g(), color.b(), h, s, v);

	h += 65.0f; // 65 instead 60 to avoid circle values
	if (h > 360.0f)
		h -= 360.0f;

	Randomizer rnd;
	s = rnd.random_float(0.65f, 1.0f);
	v = rnd.random_float(0.65f, 1.0f);

	float r, g, b;
	HSVtoRGB(h, s, v, r, g, b);
	return { r, g, b };
}

ColorRGB opposite(const ColorRGB& a, const ColorRGB& b)
{
	float ha, sa, va;
	RGBtoHSV(a.r(), a.g() , a.b(), ha, sa, va);
	float hb, sb, vb;
	RGBtoHSV(b.r(), b.g(), b.b(), hb, sb, vb);

	float delta_h = std::abs(ha - hb);
	float start_h = (delta_h > 180.0f) ? std::min(ha, hb) : std::max(ha, hb);

	start_h += 5.0f; // to avoid circle change of colors for 120 deg
	if (delta_h < 180.0f)
		delta_h = 360.0f - delta_h;

	Randomizer rnd;
	float out_h = start_h + 0.5f * delta_h;
	if (out_h > 360.0f)
		out_h -= 360.0f;
	float out_s = rnd.random_float(0.65f, 1.0f);
	float out_v = rnd.random_float(0.65f, 1.0f);

	float out_r, out_g, out_b;
	HSVtoRGB(out_h, out_s, out_v, out_r, out_g, out_b);
	return { out_r, out_g, out_b };
}

bool can_decode_color(const std::string& color) { return color.size() == 7 && color.front() == '#'; }

bool decode_color(const std::string& color_in, ColorRGB& color_out)
{
	auto hex_digit_to_int = [](const char c) {
		return
			(c >= '0' && c <= '9') ? int(c - '0') :
			(c >= 'A' && c <= 'F') ? int(c - 'A') + 10 :
			(c >= 'a' && c <= 'f') ? int(c - 'a') + 10 : -1;
	};

	color_out = ColorRGB::BLACK();
	if (can_decode_color(color_in)) {
		const char* c = color_in.data() + 1;
		for (unsigned int i = 0; i < 3; ++i) {
			const int digit1 = hex_digit_to_int(*c++);
			const int digit2 = hex_digit_to_int(*c++);
			if (digit1 != -1 && digit2 != -1)
				color_out.set(i, float(digit1 * 16 + digit2) * INV_255);
		}
	}
	else
		return false;

	assert(0.0f <= color_out.r() && color_out.r() <= 1.0f);
	assert(0.0f <= color_out.g() && color_out.g() <= 1.0f);
	assert(0.0f <= color_out.b() && color_out.b() <= 1.0f);
	return true;
}

bool decode_color(const std::string& color_in, ColorRGBA& color_out)
{
	ColorRGB rgb;
	if (!decode_color(color_in, rgb))
		return false;

	color_out = to_rgba(rgb, color_out.a());
	return true;
}

bool decode_colors(const std::vector<std::string>& colors_in, std::vector<ColorRGB>& colors_out)
{
	colors_out = std::vector<ColorRGB>(colors_in.size(), ColorRGB::BLACK());
	for (size_t i = 0; i < colors_in.size(); ++i) {
		if (!decode_color(colors_in[i], colors_out[i]))
			return false;
	}
	return true;
}

bool decode_colors(const std::vector<std::string>& colors_in, std::vector<ColorRGBA>& colors_out)
{
	colors_out = std::vector<ColorRGBA>(colors_in.size(), ColorRGBA::BLACK());
	for (size_t i = 0; i < colors_in.size(); ++i) {
		if (!decode_color(colors_in[i], colors_out[i]))
			return false;
	}
	return true;
}

std::string encode_color(const ColorRGB& color)
{
	char buffer[64];
	::sprintf(buffer, "#%02X%02X%02X", color.r_uchar(), color.g_uchar(), color.b_uchar());
	return std::string(buffer);
}

std::string encode_color(const ColorRGBA& color) { return encode_color(to_rgb(color)); }
 
ColorRGB to_rgb(const ColorRGBA& other_rgba) { return { other_rgba.r(), other_rgba.g(), other_rgba.b() }; }
ColorRGBA to_rgba(const ColorRGB& other_rgb) { return { other_rgb.r(), other_rgb.g(), other_rgb.b(), 1.0f }; }
ColorRGBA to_rgba(const ColorRGB& other_rgb, float alpha) { return { other_rgb.r(), other_rgb.g(), other_rgb.b(), alpha }; }

ColorRGBA picking_decode(unsigned int id)
{
	return {
			   float((id >> 0) & 0xff) * INV_255,  // red
			   float((id >> 8) & 0xff) * INV_255,  // green
			   float((id >> 16) & 0xff) * INV_255, // blue
			   float(picking_checksum_alpha_channel(id & 0xff, (id >> 8) & 0xff, (id >> 16) & 0xff)) * INV_255 // checksum for validating against unwanted alpha blending and multi sampling
		   };
}

unsigned int picking_encode(unsigned char r, unsigned char g, unsigned char b) { return r + (g << 8) + (b << 16); }

unsigned char picking_checksum_alpha_channel(unsigned char red, unsigned char green, unsigned char blue)
{
	// 8 bit hash for the color
	unsigned char b = ((((37 * red) + green) & 0x0ff) * 37 + blue) & 0x0ff;
	// Increase enthropy by a bit reversal
	b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
	b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
	b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
	// Flip every second bit to increase the enthropy even more.
	b ^= 0x55;
	return b;
}

Vec3f get_printable_color(std::map<int, NonplanarFacet> &triangles, Vec3f color)
{
	float mindist = 100000000.;
	std::pair<float,Vec3f> result;
	Vec3f new_color;
	//std::cout << "Ray point X:" << triangles[1].vertex->x()<< " Y:" << triangles[1].vertex->y()<< " Z:" << triangles[1].vertex->z()<< "\n";
	// find the closest point on the closest triangle of the borders
	int i = 0;
	int triangle_id = -1;
	Vec3f total_midpoint = Vec3f(0.,0.,0.);
	for(auto &facet:triangles){
		result = facet.second.find_closest_point(color);
		if(result.first < mindist){
			//std::cout << "Ray point dist:"<< result.first << " X:" << result.second.x()<< " Y:" << result.second.y()<< " Z:" << result.second.z()<< "\n";
			mindist = result.first;
			new_color = result.second;
			triangle_id = i;
		}
		// calc the triangles total midpoint
		total_midpoint += (facet.second.vertex[0]+facet.second.vertex[1]+facet.second.vertex[2])/3.;
		i++;
	}
	if(triangle_id == -1) return color;
	
	total_midpoint /= triangles.size();

	// calc each triangles midpoint
	//Vec3f triangle_midpoint = (triangles[triangle_id].vertex[0]+triangles[triangle_id].vertex[1]+triangles[triangle_id].vertex[2])/3.;
	
	// return the new point/color if the dist to total midpoint is bigger than total midpoint to triangle midpoint
	return (total_midpoint - color).norm()>(total_midpoint - new_color).norm()*1.3?new_color:color;
}

std::string ColorRGB::check_mix_color(std::vector<std::string> mix_colors)
{
	const ColorRGB tgt_color = *this;
	Vec3f tgt = Vec3f(tgt_color.r(),tgt_color.g(),tgt_color.b());
	
	std::vector<ColorRGB> mixes;
	decode_colors(mix_colors, mixes);
	// check if the color is not exactly equal to a mixer color
	for(ColorRGB & mix: mixes)
		if(mix == tgt_color) return encode_color(mix);
	//std::cout << "Color point X:" << mixes[0].r()<< " Y:" << mixes[0].g()<< " Z:" << mixes[0].b()<< "\n";
	std::map<int, NonplanarFacet> triangles;
	std::vector<int> indices;
	int n = mixes.size();
	std::vector<bool> v(n);
    for (int i = 0; i < n; ++i) {
        v[i] = (i >= (n - 3));
    }

    do {
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
               indices.push_back(i);
            }
        }
    } while (std::next_permutation(v.begin(), v.end()));
	//std::cout << "Indices point X:" << std::to_string(indices[0]) << " Y:" << std::to_string(indices[1]) << " Z:" << std::to_string(indices[2]) << "\n";
	int j = 0;
	for(int i = 0;i < indices.size(); i+=3){
		NonplanarFacet new_facet;// = new NonplanarFacet();
		new_facet.vertex[0] = {mixes[indices[i]].r(),mixes[indices[i]].g(),mixes[indices[i]].b()};
		new_facet.vertex[1] = {mixes[indices[i+1]].r(),mixes[indices[i+1]].g(),mixes[indices[i+1]].b()};
		new_facet.vertex[2] = {mixes[indices[i+2]].r(),mixes[indices[i+2]].g(),mixes[indices[i+2]].b()};
		triangles.emplace(j++,new_facet);
	}

	// check if the color is inside the available color space
	tgt = get_printable_color(triangles, tgt);

	return encode_color(ColorRGB(tgt.x(),tgt.y(),tgt.z()));	
}

std::string ColorRGB::generate_mix_gcode(std::vector<std::string> mix_colors)
{
	const ColorRGB tgt_color = *this;
	Vec3f tgt = Vec3f(tgt_color.r(),tgt_color.g(),tgt_color.b());

	std::vector<ColorRGB> mixes;
	decode_colors(mix_colors, mixes);
	std::map<int, NonplanarFacet> triangles;
	std::vector<int> indices;
	int n = mixes.size();
	std::vector<bool> v(n);
	for(int i = 0; i < mixes.size();i++)
		v.push_back(i);

	
    for (int i = 0; i < n; ++i)
        v[i] = (i >= (n - 3));
 
    do {
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                indices.push_back(i);
            }
        }
    } while (std::next_permutation(v.begin(), v.end()));
	int j = 0;
	for(int i = 0;i < indices.size(); i+=3){
		NonplanarFacet new_facet;// = new NonplanarFacet();
		new_facet.vertex[0] = {mixes[indices[i]].r(),mixes[indices[i]].g(),mixes[indices[i]].b()};
		new_facet.vertex[1] = {mixes[indices[i+1]].r(),mixes[indices[i+1]].g(),mixes[indices[i+1]].b()};
		new_facet.vertex[2] = {mixes[indices[i+2]].r(),mixes[indices[i+2]].g(),mixes[indices[i+2]].b()};
		triangles.emplace(j++,new_facet);
	}
	NonplanarSurface surf(triangles);
	// calc ratios
	// ray trace each mix color and get the right ratio for each mix color
	std::vector<float> raydists;
	std::vector<float> dists;
	std::vector<float> dists_sortable;
	float mincolmaxval = 1.;
	std::vector<float> maxvals = {0.};
	for(ColorRGB& mix: mixes){
		Vec3f dir = tgt + tgt.cross(Vec3f(mix.r(),mix.g(),mix.b())) * 10.;
		std::vector<Vec3f> dirs;
		dirs.push_back(dir);
		std::vector<IntersectInfo> hits = surf.intersect_rays(tgt,dirs);
		// find the smallest dist and use store it
		if(hits.size()>0){
			std::sort(hits.begin(), hits.end(),
              [](const auto &i1, const auto &i2) {
		  		return i2.distance > 0 && i2.bIntersect && i1.distance < i2.distance || i1.distance < 0 || !i1.bIntersect;
              });
			//std::cout <<" Y:" << hits[0].distance<< "\n"; 
		 	if(hits[0].distance>=0.) raydists.push_back(hits[0].distance);
		}
		dists.push_back((tgt - Vec3f(mix.r(),mix.g(),mix.b())).norm());
		//dists_sortable.push_back((tgt - Vec3f(mix.r(),mix.g(),mix.b())).norm());
		if (std::max(mix.r(),std::max(mix.g(),mix.b())) < mincolmaxval)
			mincolmaxval = std::max(mix.r(),std::max(mix.g(),mix.b()));
		if((tgt - Vec3f(mix.r(),mix.g(),mix.b())).norm()>maxvals[maxvals.size()-1])
			maxvals.push_back((tgt - Vec3f(mix.r(),mix.g(),mix.b())).norm());
	}
	auto above_zero = [](const float in) {
		return in<=0.?0.:in;
		};
	// build the ratios
	//std::sort(dists_sortable.begin(), dists_sortable.end());
	//float maxval = dists_sortable[dists_sortable.size()-2];
	std::vector<float> ratios;
	int i = 0;
	float sumratios = 0;
	for(float dist: dists){
		ratios.push_back(std::abs(dist==0? 1. : (1.-above_zero(std::max(raydists[i],dist*(float).99))/maxvals[maxvals.size()-1])));//std::min(360. if mincolmaxval>10 else 255. ,max(dists))
		sumratios += ratios[ratios.size()-1];
		i++;
	}
	//std::cout << "Ray point X:" << tgt_color.r()<< " Y:" << tgt_color.g()<< " Y:" << tgt_color.b()<< "\n";
	//std::cout << "Ray point X:" << sumratios << " Y:" << ratios.at(0)<< " Y:" << ratios.at(1)<< " Y:" << ratios.at(2)<< " Y:" << ratios.at(3)<< "\n";
	std::string gcode = "M165 ";
	std::string toolnrs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	i = 0;
	for(float ratio: ratios){
		// level out the ratios and export to gcode
		char rat[16];
		snprintf(rat,16, "%.2f ",(ratio>0.?ratio:0.00001)/sumratios);
		gcode += toolnrs[i++] + std::string(rat);//std::min(360. if mincolmaxval>10 else 255. ,max(dists))
	}
	//std::cout << "GCode: " << gcode << "\n";
	return gcode;
}

} // namespace Slic3r

