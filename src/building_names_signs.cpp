// 3D World - Building/Company Names and Signs
// by Frank Gennari 01/11/23

#include "function_registry.h"
#include "buildings.h"
#include "city_objects.h" // for sign_t

using std::string;

// names

string gen_random_name(rand_gen_t &rgen); // from Universe_name.cpp

string choose_family_name(rand_gen_t rgen) { // Note: deep copying so as not to update the state of the rgen that was passed in
	return gen_random_name(rgen); // use a generic random name to start with
}

namespace pixel_city {
// borrowed from Pixel City, by Shamus Young
// https://github.com/skeeto/pixelcity/blob/master/Texture.cpp
	static string const prefix[] = {"i", "Green ", "Mega", "Super ", "Omni", "e", "Hyper", "Global ", "Vital", "Next ", "Pacific ", "Metro", "Unity ", "G-",
							        "Trans", "Infinity ", "Superior ", "Monolith ", "Best ", "Atlantic ", "First ", "Union ", "National "};

	static string const name[] = {"Biotic", "Info", "Data", "Solar", "Aerospace", "Motors", "Nano", "Online", "Circuits", "Energy", "Med", "Robotic", "Exports", "Security",
		                          "Systems", "Financial", "Industrial", "Media", "Materials", "Foods", "Networks", "Shipping", "Tools", "Medical", "Publishing", "Enterprises",
		                          "Audio", "Health", "Bank", "Imports", "Apparel", "Petroleum", "Studios"};

	static string const suffix[] = {"Corp", " Inc.", "Co", "World", ".Com", " USA", " Ltd.", "Net", " Tech", " Labs", " Mfg.", " UK", " Unlimited", " One", " LLC"};

	string gen_company_name(rand_gen_t rgen) {
		string const cname(name[rgen.rand() % (sizeof(name) / sizeof(string))]);
		// randomly use a prefix OR suffix, but not both. Too verbose.
		if (rgen.rand_bool()) {return prefix[rgen.rand() % (sizeof(prefix) / sizeof(string))] + cname;}
		else                  {return cname + suffix[rgen.rand() % (sizeof(suffix) / sizeof(string))];}
	}
}

string choose_business_name(rand_gen_t rgen) {
	if (rgen.rand_bool()) {return pixel_city::gen_company_name(rgen);}
	int const v(rgen.rand()%10);

	if (v == 0) { // 3 letter acronym
		string name;
		for (unsigned n = 0; n < 3; ++n) {name.push_back('A' + rand()%26);}
		return name;
	}
	string const base(gen_random_name(rgen));
	switch (v) {
	case 1: return base;
	case 2: return base + (rgen.rand_bool() ? " Co" : " Company");
	case 3: return base + " Inc";
	case 4: return base + (rgen.rand_bool() ? " Ltd" : " Corp");
	case 5: return base + " & " + gen_random_name(rgen);
	case 6: return base + ", " + gen_random_name(rgen) + ", & " + gen_random_name(rgen);
	case 7: return (rgen.rand_bool() ? (rgen.rand_bool() ? "National " : "Global ") : (rgen.rand_bool() ? "United " : "American ")) + base;
	case 8: return base + (rgen.rand_bool() ? (rgen.rand_bool() ? " Bank" : " Trust") : (rgen.rand_bool() ? " Holdings" : " Industries"));
	case 9: return base + " " + gen_random_name(rgen);
	default: assert(0);
	}
	return ""; // never gets here
}

// signs

void gen_text_verts(vector<vert_tc_t> &verts, point const &pos, string const &text, float tsize, vector3d const &column_dir, vector3d const &line_dir, bool use_quads=0);
void add_city_building_signs(cube_t const &city_bcube, vector<sign_t> &signs);

class sign_helper_t {
	map<string, unsigned> txt_to_id;
	vector<string> text;
public:
	unsigned register_text(string const &t) {
		auto it(txt_to_id.find(t));
		if (it != txt_to_id.end()) return it->second; // found
		unsigned const id(text.size());
		txt_to_id[t] = id; // new text, insert it
		text.push_back(t);
		assert(text.size() == txt_to_id.size());
		return id;
	}
	string const &get_text(unsigned id) const {
		assert(id < text.size());
		return text[id];
	}
};

sign_helper_t sign_helper;

colorRGBA choose_sign_color(rand_gen_t &rgen) {
	colorRGBA const sign_colors[8] = {DK_RED, DK_BLUE, DK_BLUE, DK_BROWN, DK_GRAY, BLACK, BLACK, BLACK};
	return sign_colors[rgen.rand() % 8];
}

unsigned register_sign_text(string const &text) {return sign_helper.register_text(text);}

template<typename T> void add_sign_text_verts(string const &text, cube_t const &sign, bool dim, bool dir, colorRGBA const &color, vector<T> &verts_out) {
	assert(!text.empty());
	cube_t ct(sign); // text area is slightly smaller than full cube
	ct.expand_in_dim(!dim, -0.1*ct.get_sz_dim(!dim));
	ct.expand_in_dim(2, -0.05*ct.dz());
	vector3d col_dir(zero_vector), normal(zero_vector);
	bool const ldir(dim ^ dir);
	col_dir[!dim] = (ldir  ? 1.0 : -1.0);
	normal [ dim] = (dir ? 1.0 : -1.0);
	static vector<vert_tc_t> verts;
	verts.clear();
	point pos;
	pos[dim] = ct.d[dim][dir] + (dir ? 1.0 : -1.0)*0.2*ct.get_sz_dim(dim); // shift away from the front face to prevent z-fighting
	gen_text_verts(verts, pos, text, 1.0, col_dir, plus_z, 1); // use_quads=1
	assert(!verts.empty());
	cube_t text_bcube(verts[0].v);
	for (auto i = verts.begin()+2; i != verts.end(); i += 2) {text_bcube.union_with_pt(i->v);} // only need to include opposite corners
	float const width_scale(ct.get_sz_dim(!dim)/text_bcube.get_sz_dim(!dim)), height_scale(ct.dz()/text_bcube.dz());
	if (dot_product(normal, cross_product((verts[1].v - verts[0].v), (verts[2].v - verts[1].v))) < 0.0) {std::reverse(verts.begin(), verts.end());} // swap vertex winding order
	color_wrapper const cw(color);
	typename T::normal_type const nc(normal); // vector3d or norm_comp

	for (auto i = verts.begin(); i != verts.end(); ++i) {
		i->v[!dim] = i->v[!dim]*width_scale + ct.d[!dim][!ldir]; // line
		i->v.z     = i->v.z*height_scale + ct.z1(); // column
		verts_out.emplace_back(i->v, nc, i->t[0], i->t[1], cw);
	}
}
template void add_sign_text_verts(string const &text, cube_t const &sign, bool dim, bool dir, colorRGBA const &color, vector<vert_norm_comp_tc_color> &verts_out);
template void add_sign_text_verts(string const &text, cube_t const &sign, bool dim, bool dir, colorRGBA const &color, vector<vert_norm_tc_color     > &verts_out);

void add_room_obj_sign_text_verts(room_object_t const &c, colorRGBA const &color, vector<vert_norm_comp_tc_color> &verts_out) {
	add_sign_text_verts(sign_helper.get_text(c.obj_id), c, c.dim, c.dir, color, verts_out);
}

// Note: this is messy because city bcube comes from city_gen.cpp, but city buildings come from gen_buildings.cpp, neither of which has the struct definition of sign_t;
// so we need several layers of function wrappers to gather together all the pieces of data needed to add signs in this file, which includes the various required headers
void add_signs_for_city(unsigned city_id, vector<sign_t> &signs) {
	add_city_building_signs(get_city_bcube(city_id), signs);
}

void building_t::add_signs(vector<sign_t> &signs) const {
	// house welcome and other door signs are currently part of the interior - should they be? I guess at least for secondary buildings, which aren't in a city
	if (is_house) return; // no sign, for now
	if (name.empty()) return; // no company name; shouldn't get here
	rand_gen_t rgen;
	rgen.set_state((mat_ix + parts.size()*12345 + 1), (name[0] - 'A' + 1));
	rgen.rand_mix();
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // TODO: face front
	float const width(bcube.get_sz_dim(!dim)), sign_hwidth(0.25*width), sign_height(0.5*sign_hwidth), sign_depth(0.025*sign_height);
	cube_t sign;
	sign.z1() = bcube.z2(); // FIXME: top of uppermost roof
	sign.z2() = sign.z1() + sign_height;
	set_wall_width(sign, bcube.get_center_dim(!dim), sign_hwidth, !dim);
	sign.d[dim][!dir] = bcube.d[dim][ dir];
	sign.d[dim][ dir] = sign .d[dim][!dir] + (dir ? 1.0 : -1.0)*sign_depth;
	bool const two_sided(1), emissive(0);
	signs.emplace_back(sign, dim, dir, name, WHITE, choose_sign_color(rgen), two_sided, emissive);
}


