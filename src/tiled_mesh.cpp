// 3D World - Tiled Landscape Mesh Generation/Drawing Code
// by Frank Gennari
// 9/26/10

#include "GL/glew.h" // must be included first
#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"


bool const DEBUG_TILES       = 0;
int  const DISABLE_TEXTURES  = 0;
int  const TILE_RADIUS       = 4; // WM0, in mesh sizes
int  const TILE_RADIUS_IT    = 5; // WM3, in mesh sizes


extern bool tiled_mesh_display;
extern int xoff, yoff, island, DISABLE_WATER, display_mode;
extern float zmax, zmin, water_plane_z, mesh_scale;
extern point sun_pos, moon_pos;
extern float h_dirt[];
extern texture textures[];



inline float get_water_atten_factor(float mh) {
	return 2.5*WATER_COL_ATTEN*(water_plane_z - mh)*mesh_scale;
}


int get_tile_radius() {
	return ((world_mode == WMODE_INF_TERRAIN) ? TILE_RADIUS_IT : TILE_RADIUS);
}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator<(tile_xy_pair const &t) const {return ((y == t.y) ? (x < t.x) : (y < t.y));}
};


class tile_t;
tile_t *get_tile_from_xy(tile_xy_pair const &tp);


class tile_t {

	int x1, y1, x2, y2, init_dxoff, init_dyoff;
	unsigned tid, vbo, ivbo, size, stride, zvsize, base_tsize, gen_tsize;
	float radius, mzmin, mzmax, xstart, ystart, xstep, ystep;
	vector<float> zvals;
	vector<unsigned char> smask[NUM_LIGHT_SRC];
	vector<float> sh_out[NUM_LIGHT_SRC][2];

public:
	tile_t() : tid(0), vbo(0), ivbo(0), size(0), stride(0), zvsize(0), gen_tsize(0) {}
	~tile_t() {clear_vbo_tid(1,1);}
	
	tile_t(unsigned size_, int x, int y) : init_dxoff(xoff - xoff2), init_dyoff(yoff - yoff2),
		tid(0), vbo(0), ivbo(0), size(size_), stride(size+1), zvsize(stride+1), gen_tsize(0)
	{
		assert(size > 0);
		x1 = x*size;
		y1 = y*size;
		x2 = x1 + size;
		y2 = y1 + size;
		calc_start_step(0, 0);
		radius = 0.5*sqrt(xstep*xstep + ystep*ystep)*size; // approximate (lower bound)
		mzmin  = mzmax = get_camera_pos().z;
		base_tsize = get_norm_texels();
		
		if (DEBUG_TILES) {
			cout << "create " << size << ": " << x << "," << y << ", coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
		}
	}
	float get_zmin() const {return mzmin;}
	float get_zmax() const {return mzmax;}
	
	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5*(mzmin + mzmax));
	}

	void calc_start_step(int dx, int dy) {
		xstart = get_xval(x1 + dx);
		ystart = get_yval(y1 + dy);
		xstep  = (get_xval(x2 + dx) - xstart)/size;
		ystep  = (get_yval(y2 + dy) - ystart)/size;
	}

	unsigned get_gpu_memory() const {
		unsigned mem(0);
		if (vbo  > 0) mem += 2*stride*size*sizeof(vert_norm); // 44MB
		if (ivbo > 0) mem += size*size*sizeof(unsigned short); // 1MB
		if (tid  > 0) mem += 3*gen_tsize*gen_tsize; // 34MB
		return mem;
	}

	void clear_shadows() {
		for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) {
			smask[l].clear();
			for (unsigned d = 0; d < 2; ++d) sh_out[l][d].clear();
		}
	}

	void clear_tid() {
		free_texture(tid);
		gen_tsize = 0;
	}

	void clear_vbo_tid(bool vclear, bool tclear) {
		if (vclear) {
			delete_vbo(vbo);
			delete_vbo(ivbo);
			vbo = ivbo = 0;
			clear_shadows();
		}
		if (tclear) {clear_tid();}
	}

	void create_zvals() {
		RESET_TIME;
		zvals.resize(zvsize*zvsize);
		static vector<float> xv, yv; // move somewhere else?
		xv.resize(zvsize);
		yv.resize(zvsize);
		calc_start_step(0, 0);
		mzmin =  FAR_CLIP;
		mzmax = -FAR_CLIP;

		for (unsigned i = 0; i < zvsize; ++i) { // not sure about this -x, +y thing
			xv[i] = (xstart + (i - 0.5)*xstep);
			yv[i] = (ystart + (i + 0.5)*ystep);
		}
		build_xy_mesh_arrays(&xv.front(), &yv.front(), zvsize, zvsize);

		for (unsigned y = 0; y < zvsize; ++y) {
			for (unsigned x = 0; x < zvsize; ++x) {
				float const zval(fast_eval_from_index(x, y, 0, 1));
				zvals[y*zvsize + x] = zval;
				mzmin = min(mzmin, zval);
				mzmax = max(mzmax, zval);
			}
		}
		assert(mzmin <= mzmax);
		radius = 0.5*sqrt((xstep*xstep + ystep*ystep)*size*size + (mzmax - mzmin)*(mzmax - mzmin));
		//PRINT_TIME("Create Zvals");
	}

	inline vector3d get_norm(unsigned ix) const {
		return vector3d(DY_VAL*(zvals[ix] - zvals[ix + 1]), DX_VAL*(zvals[ix] - zvals[ix + zvsize]), dxdy).get_norm();
	}

	void calc_shadows(bool calc_sun, bool calc_moon) {
		bool calc_light[NUM_LIGHT_SRC] = {0};
		calc_light[LIGHT_SUN ] = calc_sun;
		calc_light[LIGHT_MOON] = calc_moon;
		tile_xy_pair const tp(x1/(int)size, y1/(int)size);

		for (unsigned l = 0; l < NUM_LIGHT_SRC; ++l) { // calculate mesh shadows for each light source
			if (!calc_light[l])    continue; // light not enabled
			if (!smask[l].empty()) continue; // already calculated (cached)
			float const *sh_in[2] = {0, 0};
			point const lpos(get_light_pos(l));
			tile_xy_pair const adj_tp[2] = {tile_xy_pair((tp.x + ((lpos.x < 0.0) ? -1 : 1)), tp.y),
				                            tile_xy_pair(tp.x, (tp.y + ((lpos.y < 0.0) ? -1 : 1)))};

			for (unsigned d = 0; d < 2; ++d) { // d = tile adjacency dimension, shared edge is in !d
				sh_out[l][!d].resize(zvsize, MESH_MIN_Z); // init value really should not be used, but it sometimes is
				tile_t *adj_tile(get_tile_from_xy(adj_tp[d]));
				if (adj_tile == NULL) continue; // no adjacent tile
				vector<float> const &adj_sh_out(adj_tile->sh_out[l][!d]);
				
				if (adj_sh_out.empty()) { // adjacent tile not initialized
					adj_tile->calc_shadows((l == LIGHT_SUN), (l == LIGHT_MOON)); // recursive call on adjacent tile
				}
				assert(adj_sh_out.size() == zvsize);
				sh_in[!d] = &adj_sh_out.front(); // chain our input to our neighbor's output
			}
			smask[l].resize(zvals.size());
			calc_mesh_shadows(l, lpos, &zvals.front(), &smask[l].front(), zvsize, zvsize,
				sh_in[0], sh_in[1], &sh_out[l][0].front(), &sh_out[l][1].front());
		}
	}

	void create_data(vector<vert_norm> &data, vector<unsigned short> &indices) {
		RESET_TIME;
		assert(zvals.size() == zvsize*zvsize);
		data.resize(stride*stride);
		indices.resize(4*size*size);
		calc_start_step(init_dxoff, init_dyoff);
		bool const has_sun(light_factor >= 0.4), has_moon(light_factor <= 0.6);
		assert(has_sun || has_moon);
		calc_shadows(has_sun, has_moon);
		
		for (unsigned y = 0; y <= size; ++y) {
			for (unsigned x = 0; x <= size; ++x) {
				unsigned const ix(y*zvsize + x);
				float light_scale(1.0);

				if (has_sun && has_moon) {
					bool const no_sun( (smask[LIGHT_SUN ][ix] & SHADOWED_ALL) != 0);
					bool const no_moon((smask[LIGHT_MOON][ix] & SHADOWED_ALL) != 0);
					light_scale = blend_light(light_factor, !no_sun, !no_moon);
				}
				else if (smask[has_sun ? LIGHT_SUN : LIGHT_MOON][ix] & SHADOWED_ALL) {
					light_scale = 0.0;
				}
				point const v((xstart + x*xstep), (ystart + y*ystep), zvals[ix]);
				data[y*stride + x].assign(v, get_norm(ix)*light_scale);
			}
		}
		for (unsigned y = 0; y < size; ++y) {
			for (unsigned x = 0; x < size; ++x) {
				unsigned const vix(y*stride + x), iix(4*(y*size + x));
				indices[iix+0] = vix;
				indices[iix+1] = vix + stride;
				indices[iix+2] = vix + stride + 1;
				indices[iix+3] = vix + 1;
			}
		}
		//PRINT_TIME("Create Data");
	}

//#define GET_TEX_DATA(t, td, id, tx, ty, tex_bs) (t[id].data + t[id].ncolors*(((ty<<tex_bs)&(t[id].height-1))*t[id].width + ((tx<<tex_bs)&(t[id].width-1))))
#define GET_TEX_DATA(t, td, id, tx, ty, tex_bs) (td[id] + t[id].ncolors*((ty&((t[id].height>>tex_bs)-1))*(t[id].width>>tex_bs) + (tx&((t[id].width>>tex_bs)-1))))

	void create_texture(unsigned tex_bs) {
		assert(tid == 0);
		assert(!island);
		assert(zvals.size() == zvsize*zvsize);
		RESET_TIME;
		unsigned const tsize(base_tsize >> tex_bs), scale(tsize/size);
		float const dz_inv(1.0/(zmax - zmin)), fscale_inv(1.0/scale);
		unsigned char *data(new unsigned char[3*tsize*tsize]); // RGB
		int k1, k2, k3, k4;
		float t;
		unsigned char const *tex_data[NUM_TEXTURES] = {0};

		for (unsigned i = 0; i < NTEX_DIRT; ++i) {
			assert(lttex_dirt[i].id < NUM_TEXTURES);
			tex_data[lttex_dirt[i].id] = textures[lttex_dirt[i].id].get_mipmap_data(tex_bs);
		}
		for (unsigned y = 0; y < size; ++y) { // makes a big performance improvement
			for (unsigned x = 0; x < size; ++x) {
				unsigned const ix(y*zvsize + x);
				float const vnz00(get_norm(ix).z), vnz01(get_norm(ix+1).z);
				float const vnz10(get_norm(ix+zvsize).z), vnz11(get_norm(ix+zvsize+1).z);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const dist(fabs(mh01 - mh00) + fabs(mh00 - mh10) + fabs(mh01 - mh11) + fabs(mh10 - mh11));
				float const relh1((min(min(mh00, mh01), min(mh10, mh11)) - zmin)*dz_inv);
				float const relh2((max(max(mh00, mh01), max(mh10, mh11)) - zmin)*dz_inv);
				get_tids(relh1, NTEX_DIRT-1, h_dirt, k1, k2, t);
				get_tids(relh2, NTEX_DIRT-1, h_dirt, k3, k4, t);
				bool const same_tid(k1 == k4);
				k2 = k4;
				
				for (unsigned sy = 0; sy < scale; ++sy) {
					unsigned const ty(y*scale + sy);
					float const ypi(sy*fscale_inv);

					for (unsigned sx = 0; sx < scale; ++sx) {
						unsigned const tx(x*scale + sx), ix(y*zvsize + x), off(3*(ty*tsize + tx));
						float const xpi(sx*fscale_inv);
						float const mh((1.0 - xpi)*((1.0 - ypi)*mh00 + ypi*mh10) + xpi*((1.0 - ypi)*mh01 + ypi*mh11));

						if (!same_tid) {
							float const relh((mh - zmin)*dz_inv);
							get_tids(relh, NTEX_DIRT-1, h_dirt, k1, k2, t);
						}
						int const id(lttex_dirt[k1].id), id2(lttex_dirt[k2].id);
						unsigned char const *t1_data(GET_TEX_DATA(textures, tex_data, id, tx, ty, tex_bs));
						unsigned char *td(data + off);
						
						if (k1 == k2) { // single texture
							RGB_BLOCK_COPY(td, t1_data);
						}
						else { // blend two textures - performance critical
							unsigned char const *t2_data(GET_TEX_DATA(textures, tex_data, id2, tx, ty, tex_bs));
							BLEND_COLOR(td, t2_data, t1_data, t);
						}

						// handle steep slopes (dirt/rock texture replaces grass texture)
						float const sthresh[2] = {0.45, 0.7};
						float const vnz((1.0 - xpi)*((1.0 - ypi)*vnz00 + ypi*vnz10) + xpi*((1.0 - ypi)*vnz01 + ypi*vnz11));

						if (vnz < sthresh[1]) {
							if (id == GROUND_TEX || id2 == GROUND_TEX) { // ground/grass
								assert(tex_data[DIRT_TEX]);
								unsigned char const *ta_data(GET_TEX_DATA(textures, tex_data, DIRT_TEX, tx, ty, tex_bs));
								unsigned char temp[3];

								if (id == GROUND_TEX || id2 == ROCK_TEX) {
									assert(tex_data[ROCK_TEX]);
									unsigned char const *tb_data(GET_TEX_DATA(textures, tex_data, ROCK_TEX, tx, ty, tex_bs));
									BLEND_COLOR(temp, tb_data, ta_data, t);
								}
								else {
									RGB_BLOCK_COPY(temp, ta_data);
								}
								float const val(CLIP_TO_01((vnz - sthresh[0])/(sthresh[1] - sthresh[0])));
								BLEND_COLOR(td, td, temp, val);
							}
							else if (id2 == SNOW_TEX) { // snow
								assert(tex_data[ROCK_TEX]);
								unsigned char const *ta_data(GET_TEX_DATA(textures, tex_data, ROCK_TEX, tx, ty, tex_bs));
								float const val(CLIP_TO_01(2.0f*(vnz - sthresh[0])/(sthresh[1] - sthresh[0])));
								BLEND_COLOR(td, td, ta_data, val);
							}
						}

						// darken underwater regions
						if (!DISABLE_WATER && mh < water_plane_z) {
							float c[3] = {td[0], td[1], td[2]};
							atten_by_water_depth(c, get_water_atten_factor(mh));
							UNROLL_3X(td[i_] = (unsigned char)c[i_];)
						}
					} // for sx
				} // for sy
			} // for x
		} // for y
		//if (tex_bs == 0) {PRINT_TIME("Texture Gen");}
		bool const mipmaps(0); // mipmaps take about 33% more memory
		setup_texture(tid, GL_MODULATE, mipmaps, 0, 0, 0, 0);
		assert(tid > 0);
		assert(glIsTexture(tid));
		//bool const has_comp(has_extension("GL_ARB_texture_compression")); GL_COMPRESSED_RGB - too slow
		//if (mipmaps) glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, tsize, tsize, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		if (mipmaps) gen_mipmaps();
		glDisable(GL_TEXTURE_2D);
		delete [] data;
		//if (tex_bs == 0) {PRINT_TIME("Texture Upload");}
	}

	void check_texture() {
		float dist(get_rel_dist_to_camera()*get_tile_radius()); // in tiles
		unsigned tex_bs(0);

		while (dist > 1.0 && (base_tsize >> tex_bs) > size) { // texture LOD
			dist /= 2;
			++tex_bs;
		}
		unsigned const new_tsize(base_tsize >> tex_bs);

		//if (new_tsize > gen_tsize) {
		//if (new_tsize != gen_tsize) { // power of 2 >= size
		if (new_tsize > gen_tsize || 2*new_tsize < gen_tsize) {
			clear_tid();
			gen_tsize = new_tsize; // power of 2 >= size
			create_texture(tex_bs);
		}
	}

	float get_rel_dist_to_camera() const {
		return (p2p_dist_xy(get_camera_pos(), get_center()) - radius)/(get_tile_radius()*(X_SCENE_SIZE + Y_SCENE_SIZE));
	}

	bool update_range() { // if returns 0, tile will be deleted
		float const dist(get_rel_dist_to_camera());
		if (dist > 1.5) clear_vbo_tid(1,1);
		return (dist < 2.0);
	}

	bool is_visible() const {
		return camera_pdu.sphere_visible_test(get_center(), radius);
	}

	void bind_vbos() {
		assert(vbo > 0 && ivbo > 0);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo, 1);
	}

	bool draw(vector<vert_norm> &data, vector<unsigned short> &indices) { // make const or make vbo mutable?
		assert(size > 0);
		if (!is_visible()) return 0; // not visible to camera

		if (!DISABLE_TEXTURES) {
			check_texture();
			
			if (tid > 0) {
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, tid);
			}
		}
		glPushMatrix();
		glTranslatef(((xoff - xoff2) - init_dxoff)*DX_VAL, ((yoff - yoff2) - init_dyoff)*DY_VAL, 0.0);
		if (tid > 0) set_landscape_texgen(1.0, (-x1 - init_dxoff), (-y1 - init_dyoff), MESH_X_SIZE, MESH_Y_SIZE);
		unsigned ptr_stride(sizeof(vert_norm));

		if (vbo == 0) {
			create_data(data, indices);
			vbo  = create_vbo();
			ivbo = create_vbo();
			bind_vbos();
			upload_vbo_data(&data.front(),    data.size()*ptr_stride,                0);
			upload_vbo_data(&indices.front(), indices.size()*sizeof(unsigned short), 1);
		}
		else {
			bind_vbos();
		}
		glVertexPointer(3, GL_FLOAT, ptr_stride, 0);
		glNormalPointer(   GL_FLOAT, ptr_stride, (void *)sizeof(point));
		glDrawRangeElements(GL_QUADS, 0, data.size(), 4*size*size, GL_UNSIGNED_SHORT, 0); // requires GL/glew.h
		//glDrawElements(GL_QUADS, 4*size*size, GL_UNSIGNED_SHORT, 0);
		bind_vbo(0, 0);
		bind_vbo(0, 1);
		glPopMatrix();
		if (tid > 0) disable_textures_texgen();
		return 1;
	}
};


class tile_draw_t {

	typedef map<tile_xy_pair, tile_t*> tile_map;
	tile_map tiles;

public:
	tile_draw_t() {
		assert(MESH_X_SIZE == MESH_Y_SIZE && X_SCENE_SIZE == Y_SCENE_SIZE);
	}
	~tile_draw_t() {clear();}

	void clear() {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			delete i->second;
		}
		tiles.clear();
	}

	void update() {
		//RESET_TIME;
		assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
		point const camera(get_camera_pos() - point((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
		int const tile_radius(int(1.5*get_tile_radius()) + 1);
		int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
		int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
		int const x2( tile_radius + toffx), y2( tile_radius + toffy);
		unsigned const init_tiles(tiles.size());
		vector<tile_xy_pair> to_erase;

		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) { // free old tiles
			if (!i->second->update_range()) {
				to_erase.push_back(i->first);
				delete i->second;
			}
		}
		for (vector<tile_xy_pair>::const_iterator i = to_erase.begin(); i != to_erase.end(); ++i) {
			tiles.erase(*i);
		}
		for (int y = y1; y <= y2; ++y ) { // create new tiles
			for (int x = x1; x <= x2; ++x ) {
				tile_xy_pair const txy(x, y);
				if (tiles.find(txy) != tiles.end()) continue; // already exists
				tile_t tile(MESH_X_SIZE, x, y);
				if (tile.get_rel_dist_to_camera() >= 1.5) continue; // too far away to create
				tile_t *new_tile(new tile_t(tile));
				new_tile->create_zvals();
				tiles[txy] = new_tile;
			}
		}
		if (DEBUG_TILES && (tiles.size() != init_tiles || !to_erase.empty())) {
			cout << "update: tiles: " << init_tiles << " to " << tiles.size() << ", erased: " << to_erase.size() << endl;
		}
		//PRINT_TIME("Tiled Terrain Update");
	}

	float draw(bool add_hole) {
		float zmin(FAR_CLIP);
		glDisable(GL_NORMALIZE);
		glDisable(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisable(GL_COLOR_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		setup_mesh_lighting();
		unsigned num_drawn(0);
		unsigned long long mem(0);
		vector<vert_norm> data;
		vector<unsigned short> indices;
		static point last_sun(all_zeros), last_moon(all_zeros);
		static bool last_water_en(1);
		bool const water_en((display_mode & 0x04) != 0);

		if (sun_pos != last_sun || moon_pos != last_moon) {
			clear_vbos_tids(1,0); // light source changed, clear vbos and build new shadow map
			last_sun  = sun_pos;
			last_moon = moon_pos;
		}
		if (water_en != last_water_en) { // should this be here?
			clear_vbos_tids(0,1); // water changed, recreate textures
			last_water_en = water_en;
		}
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			assert(i->second);
			if (DEBUG_TILES) mem += i->second->get_gpu_memory();
			if (add_hole && i->first.x == 0 && i->first.y == 0) continue;
			if (i->second->get_rel_dist_to_camera() > 1.4)      continue; // too far to draw
			zmin = min(zmin, i->second->get_zmin());
			num_drawn += i->second->draw(data, indices);
		}
		if (DEBUG_TILES) cout << "tiles drawn: " << num_drawn << " of " << tiles.size() << ", gpu mem: " << mem/1024/1024 << endl;
		run_post_mesh_draw();
		return zmin;
	}

	void clear_vbos_tids(bool vclear, bool tclear) {
		for (tile_map::iterator i = tiles.begin(); i != tiles.end(); ++i) {
			i->second->clear_vbo_tid(vclear, tclear);
		}
	}

	tile_t *get_tile_from_xy(tile_xy_pair const &tp) {
		tile_map::iterator it(tiles.find(tp));
		if (it != tiles.end()) return it->second;
		return NULL;
	}
};


tile_draw_t terrain_tile_draw;


tile_t *get_tile_from_xy(tile_xy_pair const &tp) {

	return terrain_tile_draw.get_tile_from_xy(tp);
}


void draw_vert_color(colorRGBA c, float x, float y, float z) {

	if (z < water_plane_z) atten_by_water_depth(&c.red, get_water_atten_factor(z));
	c.do_glColor();
	glVertex3f(x, y, z);
}


void fill_gap() {

	//RESET_TIME;
	colorRGBA const color(setup_mesh_lighting());
	
	if (!DISABLE_TEXTURES) {
		select_texture(LANDSCAPE_TEX);
		set_landscape_texgen(1.0, xoff, yoff, MESH_X_SIZE, MESH_Y_SIZE);
	}
	vector<float> xv(MESH_X_SIZE+1), yv(MESH_Y_SIZE+1);
	float const xstart(get_xval(xoff2)), ystart(get_yval(yoff2));

	for (int i = 0; i <= MESH_X_SIZE; ++i) { // not sure about this -x, +y thing
		xv[i] = (xstart + (i - 0.5)*DX_VAL);
	}
	for (int i = 0; i <= MESH_Y_SIZE; ++i) {
		yv[i] = (ystart + (i + 0.5)*DY_VAL);
	}

	// draw +x
	build_xy_mesh_arrays(&xv.front(), &yv[MESH_Y_SIZE], MESH_X_SIZE, 1);
	glBegin(GL_QUAD_STRIP);

	for (int x = 0; x < MESH_X_SIZE; ++x) {
		vertex_normals[MESH_Y_SIZE-1][x].do_glNormal();
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE-DY_VAL, mesh_height[MESH_Y_SIZE-1][x]);
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE       , fast_eval_from_index(x, 0, 0, 1));
	}
	glEnd();
	float const z_end_val(fast_eval_from_index(MESH_X_SIZE-1, 0, 0, 1));

	// draw +y
	build_xy_mesh_arrays(&xv[MESH_X_SIZE], &yv.front(), 1, MESH_Y_SIZE+1);
	glBegin(GL_QUAD_STRIP);

	for (int y = 0; y < MESH_X_SIZE; ++y) {
		vertex_normals[y][MESH_X_SIZE-1].do_glNormal();
		draw_vert_color(color, X_SCENE_SIZE-DX_VAL, get_yval(y), mesh_height[y][MESH_X_SIZE-1]);
		draw_vert_color(color, X_SCENE_SIZE       , get_yval(y), fast_eval_from_index(0, y, 0, 1));
	}

	// draw corner quad
	draw_vert_color(color, X_SCENE_SIZE-DX_VAL, Y_SCENE_SIZE, z_end_val);
	draw_vert_color(color, X_SCENE_SIZE       , Y_SCENE_SIZE, fast_eval_from_index(0, MESH_Y_SIZE, 0, 1));
	glEnd();

	disable_textures_texgen();
	glDisable(GL_COLOR_MATERIAL);
	run_post_mesh_draw();
	//PRINT_TIME("Fill Gap");
}


float draw_tiled_terrain(bool add_hole) {

	//RESET_TIME;
	bool const vbo_supported(setup_gen_buffers());
		
	if (!vbo_supported) {
		cout << "Warning: VBOs not supported, so tiled mesh cannot be enabled." << endl;
		tiled_mesh_display = 0;
		return zmin;
	}
	terrain_tile_draw.update();
	float const zmin(terrain_tile_draw.draw(add_hole));
	if (add_hole) fill_gap(); // need to fill the gap on +x/+y
	//glFinish(); PRINT_TIME("Tiled Terrain Draw"); //exit(0);
	return zmin;
}


void clear_tiled_terrain() {
	terrain_tile_draw.clear();
}

void reset_tiled_terrain_state() {
	terrain_tile_draw.clear_vbos_tids(1,1);
}


