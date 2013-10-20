// 3D World - Heightmap Texture Managment
// by Frank Gennari
// 10/19/13

#include "heightmap.h"
#include "function_registry.h"
#include "inlines.h"

using namespace std;


unsigned const TEX_EDGE_MODE = 2; // 0 = clamp, 1 = cliff/underwater, 2 = mirror

extern float mesh_scale, dxdy;

float scale_mh_texture_val(float val);


float heightmap_t::get_heightmap_value(unsigned x, unsigned y) const { // returns values from 0 to 256

	assert(is_allocated());
	assert(ncolors == 1 || ncolors == 2); // one or two byte grayscale
	assert(x < (unsigned)width && y < (unsigned)height);
	unsigned const ix(width*y + x);
	return ((ncolors == 2) ? (data[ix<<1]/256.0 + data[(ix<<1)+1]) : data[ix]);
}


void heightmap_t::modify_heightmap_value(unsigned x, unsigned y, int val, bool val_is_delta) {

	assert(is_allocated());
	assert(ncolors == 1 || ncolors == 2); // one or two byte grayscale
	assert(x < (unsigned)width && y < (unsigned)height);
	unsigned const ix(width*y + x);

	if (ncolors == 1) {
		if (val_is_delta) {val += data[ix];}
		data[ix] = max(0, min(255, val)); // clamp
	}
	else { // ncolors == 2
		unsigned short *ptr((unsigned short *)(data[ix<<1]));
		if (val_is_delta) {val += *ptr;}
		*ptr = max(0, min(65535, val)); // clamp
	}
}


void tex_mod_map_manager_t::add_mod(tex_mod_vect_t const &mod) { // vector (could use a template function)
	for (tex_mod_vect_t::const_iterator i = mod.begin(); i != mod.end(); ++i) {add_mod(*i);}
}

void tex_mod_map_manager_t::add_mod(tex_mod_map_t const &mod) { // map (could use a template function)
	for (tex_mod_map_t::const_iterator i = mod.begin(); i != mod.end(); ++i) {add_mod(mod_elem_t(*i));}
}

bool tex_mod_map_manager_t::read_mod(string const &fn) {

	//assert(mod_map.empty()); // ???
	mod_map.clear(); // allow merging ???
	FILE *fp(fopen(fn.c_str(), "rb"));

	if (fp == NULL) {
		cerr << "Error opening terrain height mod map " << fn << " for read" << endl;
		return 0;
	}
	unsigned sz(0);
	unsigned const sz_read(fread(&sz, sizeof(unsigned), 1, fp));
	assert(sz_read == 1); // add error checking?

	for (unsigned i = 0; i < sz; ++i) {
		mod_elem_t elem;
		unsigned const elem_read(fread(&elem, sizeof(mod_elem_t), 1, fp)); // use a larger block?
		assert(elem_read == 1); // add error checking?
		mod_map.add(elem);
	}
	fclose(fp);
	return 1;
}

bool tex_mod_map_manager_t::write_mod(string const &fn) const {

	FILE *fp(fopen(fn.c_str(), "wb"));

	if (fp == NULL) {
		cerr << "Error opening terrain height mod map " << fn << " for write" << endl;
		return 0;
	}
	unsigned const sz(mod_map.size());
	unsigned const sz_write(fwrite(&sz, sizeof(unsigned), 1, fp));
	assert(sz_write == 1); // add error checking?

	for (tex_mod_map_t::const_iterator i = mod_map.begin(); i != mod_map.end(); ++i) {
		mod_elem_t const elem(*i); // could use *i directly?
		unsigned const elem_write(fwrite(&elem, sizeof(mod_elem_t), 1, fp)); // use a larger block?
		assert(elem_write == 1); // add error checking?
	}
	fclose(fp);
	return 1;
}


bool terrain_hmap_manager_t::clamp_xy(int &x, int &y) const {

	assert(hmap.width > 0 && hmap.height > 0);
	x = round_fp(mesh_scale*x) + hmap.width /2; // scale and offset (0,0) to texture center
	y = round_fp(mesh_scale*y) + hmap.height/2;

	switch (TEX_EDGE_MODE) {
	case 0: // clamp
		x = max(0, min(hmap.width -1, x));
		y = max(0, min(hmap.height-1, y));
		break;
	case 1: // cliff/underwater
		if (x < 0 || y < 0 || x >= hmap.width || y >= hmap.height) {return 0;} // off the texture
		break;
	case 2: // mirror
		{
			int const xmod(abs(x)%hmap.width), ymod(abs(y)%hmap.height), xdiv(x/hmap.width), ydiv(y/hmap.height);
			x = ((xdiv & 1) ? (hmap.width  - xmod - 1) : xmod);
			y = ((ydiv & 1) ? (hmap.height - ymod - 1) : ymod);
		}
		break;
	}
	return 1;
}

bool terrain_hmap_manager_t::clamp_elem(mod_elem_t &elem) {

	assert((unsigned)max(hmap.width, hmap.height) <= max_tex_ix());
	int x(elem.x), y(elem.y);
	if (!clamp_xy(x, y)) return 0;
	assert(x >= 0 && y >= 0 && x < hmap.width && y < hmap.height);
	elem.x = x; elem.y = y;
	return 1;
}

void terrain_hmap_manager_t::load(char const *const fn, bool invert_y) {

	assert(fn != NULL);
	cout << "Loading terrain heightmap file " << fn << endl;
	RESET_TIME;
	assert(!hmap.is_allocated()); // can only call once
	hmap = heightmap_t(0, 7, 0, 0, fn, invert_y);
	hmap.load(-1, 0, 1, 1);
	PRINT_TIME("Heightmap Load");
}

float terrain_hmap_manager_t::get_clamped_height(int x, int y) const { // translate so that (0,0) is in the center of the heightmap texture

	assert(enabled());
	if (!clamp_xy(x, y)) {return scale_mh_texture_val(0.0);} // off the texture, use min value
	return scale_mh_texture_val(hmap.get_heightmap_value(x, y));
}

float terrain_hmap_manager_t::interpolate_height(float x, float y) const { // bilinear interpolation

	int const xlo(floor(x)), ylo(floor(y));
	float const xv(x - xlo), yv(y - ylo);
	return   yv *(xv*get_clamped_height(xlo+1, ylo+1) + (1.0-xv)*get_clamped_height(xlo, ylo+1)) +
		(1.0-yv)*(xv*get_clamped_height(xlo+1, ylo  ) + (1.0-xv)*get_clamped_height(xlo, ylo  ));
}

vector3d terrain_hmap_manager_t::get_norm(int x, int y) const {
	return vector3d(DY_VAL*(get_clamped_height(x, y) - get_clamped_height(x+1, y)),
			        DX_VAL*(get_clamped_height(x, y) - get_clamped_height(x, y+1)), dxdy).get_norm();
}

bool terrain_hmap_manager_t::modify_height(mod_elem_t &elem) { // elem may be modified

	if (!clamp_elem(elem)) return 0;
	hmap.modify_heightmap_value(elem.x, elem.y, elem.delta, 1);
	return 1;
}

bool terrain_hmap_manager_t::modify_and_cache_height(mod_elem_t &elem) { // elem may be modified

	if (!modify_height(elem)) return 0;
	add_mod(elem);
	return 1;
}

bool terrain_hmap_manager_t::read_mod(string const &fn) {
	if (!tex_mod_map_manager_t::read_mod(fn)) return 0;

	for (tex_mod_map_t::const_iterator i = mod_map.begin(); i != mod_map.end(); ++i) { // apply the mod to the current texture
		assert(i->first.x < hmap.width && i->first.y < hmap.height); // ensure the mod values fit within the texture
		bool const ret(modify_height(mod_elem_t(*i))); // should *not* modify the temporary
		assert(ret); // require success
	}
	return 1;
}


