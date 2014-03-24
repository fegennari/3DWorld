// 3D World - Drawing Code
// by Frank Gennari
// 3/10/02
#include "3DWorld.h"
#include "function_registry.h"
#include "shaders.h"


string const default_font_texture_atlas_fn = "textures/atlas/text_atlas.png";
string font_texture_atlas_fn(default_font_texture_atlas_fn);

extern int window_height, display_mode;


struct per_char_data_t {
	float u1, u2, v1, v2; // texture coords
	float width; // horizontal size/step size

	per_char_data_t(float u1_=0, float u2_=0, float v1_=0, float v2_=0, float width_=1) : u1(u1_), u2(u2_), v1(v1_), v2(v2_), width(width_) {}
};


class font_texture_manager_t {

	texture_t texture;
	per_char_data_t pcd[256]; // one per ASCII value

public:
	font_texture_manager_t() : texture(0, 7, 0, 0, 0, 4, 3, "", 0, 0) {} // custom alpha mipmaps, uncompressed

	bool check_nonempty_tile_column(unsigned tx, unsigned ty, unsigned x, unsigned tsize) const {
		unsigned char const *data(texture.get_data());

		for (unsigned y = 0; y < tsize; ++y) { // iterate over one character tile
			unsigned const pixel((ty*tsize + y)*texture.width + tx*tsize + x);
			if (data[(pixel<<2)+3] != 0) return 1; // check alpha component
		}
		return 0;
	}

	// Note: expects to find a square texture with 16x16 tiles, one per ASCII value, starting from the ULC
	void load(string const &fn) {
		texture.free_data();
		texture.name = (fn.empty() ? font_texture_atlas_fn : fn);
		texture.load(-1);
		assert(texture.ncolors == 4); // RGBA (really only need RA though)
		assert(texture.width == texture.height);
		assert((texture.width & 15) == 0); // multiple of 16
		float const duv(1.0/16.0), pw(1.0/texture.width);
		unsigned const tsize(texture.width >> 4); // tile size
		unsigned char *data(texture.get_data());

		for (unsigned i = 0; i < texture.num_pixels(); ++i) {
			unsigned weight(0);
			UNROLL_3X(weight += data[4*i+i_];)
			data[4*i+3] = (unsigned char)(weight/3.0); // set alpha equal to luminance
			UNROLL_3X(data[4*i+i_] = 255;) // set luminance/RGB to 1
		}
		for (unsigned ty = 0; ty < 16; ++ty) {
			for (unsigned tx = 0; tx < 16; ++tx) {
				// calculate kerning by looking for all black/transparent columns
				unsigned col_start(0), col_end(tsize/2); // non-printable chars are half width
				unsigned const ty_inv(15-ty), dix((ty_inv<<4) + tx); // index into pcd

				for (unsigned x = 0; x < tsize; ++x) { // forward iterate over one character tile
					if (check_nonempty_tile_column(tx, ty, x, tsize)) {col_start = x; break;}
				}
				for (unsigned x = tsize; x > col_start; --x) { // backward iterate over one character tile
					if (check_nonempty_tile_column(tx, ty, x-1, tsize)) {col_end = x; break;}
				}
				float const width(float(col_end - col_start)/float(tsize));
				pcd[dix] = per_char_data_t(tx*duv+pw*col_start, (tx+1)*duv-pw*(tsize - col_end), ty*duv, (ty+1)*duv-pw, width);
			} // for tx
		} // for ty
	}
	void bind_gl() {
		texture.check_init();
		texture.bind_gl();
	}
	per_char_data_t const &lookup_ascii(unsigned char val) const {assert(val < 256); return pcd[val];}
	void free_gl_state() {texture.gl_delete();}
};

font_texture_manager_t font_texture_manager; // singleton

void load_font_texture_atlas(string const &fn) {font_texture_manager.load(fn);}
void free_font_texture_atlas() {font_texture_manager.free_gl_state();}


void draw_bitmap_text(colorRGBA const &color, point const &pos, string const &text, float tsize) {

	if (text.empty()) return; // error?
	float const line_spacing = 1.25;
	float const char_spacing = 0.06;
	float const char_sz(0.001*tsize);
	point cursor(pos);
	vector<vert_tc_t> verts;

	for (string::const_iterator i = text.begin(); i != text.end(); ++i) {
		if (*i == '\n') { // newline (CR/LF)
			cursor.x  = pos.x; // CR
			cursor.y -= line_spacing*char_sz; // LF
		}
		else {
			per_char_data_t const &pcd(font_texture_manager.lookup_ascii(*i));
			if (pcd.width == 0.0) continue; // non-printable character, skip it (currently never get here, but left for future use)
			float const char_width(char_sz*pcd.width);
			float const t[4][2] = {{pcd.u1,pcd.v1}, {pcd.u2,pcd.v1}, {pcd.u2,pcd.v2}, {pcd.u1,pcd.v2}};
			float const dx[4] = {0.0, char_width, char_width, 0.0};
			float const dy[4] = {0.0, 0.0, char_sz, char_sz};

			for (unsigned i = 0; i < 6; ++i) {
				unsigned const ix(quad_to_tris_ixs[i]);
				point p(cursor);
				p.x += dx[ix]; p.y += dy[ix];
				verts.push_back(vert_tc_t(p, t[ix][0], t[ix][1]));
			}
			cursor.x += char_width + char_sz*char_spacing;
		}
	}
	shader_t s;
	s.begin_simple_textured_shader(0.1, 0, 0, &color);
	font_texture_manager.bind_gl();
	enable_blend();
	draw_verts(verts, GL_TRIANGLES);
	disable_blend();
	s.end_shader();
}


void draw_text(colorRGBA const &color, float x, float y, float z, char const *text, float tsize) {

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // always filled
	glDisable(GL_DEPTH_TEST);
	draw_bitmap_text(color, point(x, y, z), text, 0.8*tsize);
	glEnable(GL_DEPTH_TEST);
}

