// 3D World - Tiled Landscape Mesh Generation/Drawing Code
// by Frank Gennari
// 9/26/10

#include "3DWorld.h"
#include "mesh.h"
#include "textures_3dw.h"
#include "gl_ext_arb.h"
#include "shaders.h"
#include "small_tree.h"


bool const DEBUG_TILES        = 0;
bool const ENABLE_TREE_LOD    = 1; // faster but has popping artifacts
bool const ENABLE_TREE_BFC    = 0; // faster but has popping artifacts
int  const DISABLE_TEXTURES   = 0;
int  const TILE_RADIUS        = 4; // WM0, in mesh sizes
int  const TILE_RADIUS_IT     = 6; // WM3, in mesh sizes
unsigned const NUM_LODS       = 5; // > 0
float const DRAW_DIST_TILES   = 1.4;
float const CREATE_DIST_TILES = 1.5;
float const CLEAR_DIST_TILES  = 1.5;
float const DELETE_DIST_TILES = 2.0;


extern int xoff, yoff, island, DISABLE_WATER, display_mode, show_fog, tree_mode;
extern float zmax, zmin, water_plane_z, mesh_scale, mesh_scale2, mesh_scale_z, vegetation, relh_adj_tex;
extern point sun_pos, moon_pos;
extern float h_dirt[];
extern texture_t textures[];

void draw_water_edge(float zval);


int get_tile_radius() {
	return ((world_mode == WMODE_INF_TERRAIN) ? TILE_RADIUS_IT : TILE_RADIUS);
}

float get_inf_terrain_fog_dist() {
	return DRAW_DIST_TILES*get_tile_radius()*(X_SCENE_SIZE + Y_SCENE_SIZE);
}

bool is_water_enabled() {return (!DISABLE_WATER && (display_mode & 0x04) != 0);}
bool trees_enabled   () {return (world_mode == WMODE_INF_TERRAIN && (tree_mode & 2) && vegetation > 0.0);}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator<(tile_xy_pair const &t) const {return ((y == t.y) ? (x < t.x) : (y < t.y));}
};


class tile_t;
tile_t *get_tile_from_xy(tile_xy_pair const &tp);


class tile_t {

	int x1, y1, x2, y2, init_dxoff, init_dyoff, init_tree_dxoff, init_tree_dyoff;
	unsigned tid, vbo, ivbo[NUM_LODS], size, stride, zvsize, base_tsize, gen_tsize, tree_lod_level;
	float radius, mzmin, mzmax, xstart, ystart, xstep, ystep;
	vector<float> zvals;
	vector<unsigned char> smask[NUM_LIGHT_SRC];
	vector<float> sh_out[NUM_LIGHT_SRC][2];
	small_tree_group trees;

public:
	typedef vert_norm_comp vert_type_t;

	tile_t() : tid(0), vbo(0), size(0), stride(0), zvsize(0), gen_tsize(0), tree_lod_level(0) {init_vbo_ids();}
	~tile_t() {clear_vbo_tid(1,1);}
	
	tile_t(unsigned size_, int x, int y) : init_dxoff(xoff - xoff2), init_dyoff(yoff - yoff2), init_tree_dxoff(0), init_tree_dyoff(0),
		tid(0), vbo(0), size(size_), stride(size+1), zvsize(stride+1), gen_tsize(0), tree_lod_level(0)
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
		init_vbo_ids();
		
		if (DEBUG_TILES) {
			cout << "create " << size << ": " << x << "," << y << ", coords: " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
		}
	}
	float get_zmin() const {return mzmin;}
	float get_zmax() const {return mzmax;}
	bool has_water() const {return (mzmin < water_plane_z);}
	
	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5*(mzmin + mzmax));
	}
	cube_t get_cube() const {
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2));
		return cube_t(xv1, xv1+(x2-x1)*DX_VAL, yv1, yv1+(y2-y1)*DY_VAL, mzmin, mzmax);
	}
	bool contains_camera() const {
		return get_cube().contains_pt_xy(get_camera_pos());
	}

	void calc_start_step(int dx, int dy) {
		xstart = get_xval(x1 + dx);
		ystart = get_yval(y1 + dy);
		xstep  = (get_xval(x2 + dx) - xstart)/size;
		ystep  = (get_yval(y2 + dy) - ystart)/size;
	}

	unsigned get_gpu_memory() const {
		unsigned mem(0);
		if (vbo > 0) mem += 2*stride*size*sizeof(vert_type_t);
		if (tid > 0) mem += 4*size*size; // 4 bytes per texel

		for (unsigned i = 0; i < NUM_LODS; ++i) {
			if (ivbo[i] > 0) mem += (size>>i)*(size>>i)*sizeof(unsigned short);
		}
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

	void init_vbo_ids() {
		vbo = 0;
		for (unsigned i = 0; i < NUM_LODS; ++i) {ivbo[i] = 0;}
	}

	void clear_vbo_tid(bool vclear, bool tclear) {
		if (vclear) {
			delete_vbo(vbo);
			for (unsigned i = 0; i < NUM_LODS; ++i) {delete_vbo(ivbo[i]);}
			init_vbo_ids();
			clear_shadows();
			trees.vbo_manager.clear_vbo();
		}
		if (tclear) {clear_tid();}
	}

	void create_xy_arrays(unsigned xy_size, float xy_scale) {
		static vector<float> xv, yv; // move somewhere else?
		xv.resize(xy_size);
		yv.resize(xy_size);

		for (unsigned i = 0; i < xy_size; ++i) { // not sure about this -x, +y thing
			xv[i] = xy_scale*(xstart + (i - 0.5)*xstep);
			yv[i] = xy_scale*(ystart + (i + 0.5)*ystep);
		}
		build_xy_mesh_arrays(&xv.front(), &yv.front(), xy_size, xy_size);
	}

	void create_zvals() {
		//RESET_TIME;
		zvals.resize(zvsize*zvsize);
		calc_start_step(0, 0);
		create_xy_arrays(zvsize, 1.0);
		mzmin =  FAR_CLIP;
		mzmax = -FAR_CLIP;

		#pragma omp parallel for schedule(static,1)
		for (int y = 0; y < (int)zvsize; ++y) {
			for (unsigned x = 0; x < zvsize; ++x) {
				zvals[y*zvsize + x] = fast_eval_from_index(x, y);
			}
		}
		for (vector<float>::const_iterator i = zvals.begin(); i != zvals.end(); ++i) {
			mzmin = min(mzmin, *i);
			mzmax = max(mzmax, *i);
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

	void apply_tree_ao_shadows(vector<vert_type_t> &data) { // should this generate a float or unsigned char shadow weight instead?
		point const tree_off((init_dxoff - init_tree_dxoff)*DX_VAL, (init_dyoff - init_tree_dyoff)*DY_VAL, 0.0);

		for (small_tree_group::const_iterator i = trees.begin(); i != trees.end(); ++i) {
			point const pos(i->get_pos() + tree_off);
			float const radius(i->get_pine_tree_radius());
			int const xc((pos.x - xstart)/xstep), yc((pos.y - ystart)/ystep);
			if (xc < 0 || yc < 0 || xc >= (int)stride || yc >= (int)stride) continue; // off the edge of the mesh (shouldn't occur)
			int rval(max(int(radius/xstep), int(radius/ystep)) + 2);
			rval = min(rval, min(xc, yc)); // clip to lower bounds
			rval = min(rval, min((int)stride-xc-1, (int)stride-yc-1)); // clip to upper bounds
			float const rval_inv(1.0/rval);
			int const x1(xc-rval), y1(yc-rval), x2(xc+rval), y2(yc+rval);
			assert(x1 >= 0 && y1 >= 0 && x2 < (int)stride && y2 < (int)stride);

			for (int y = y1; y <= y2; ++y) {
				for (int x = x1; x <= x2; ++x) {
					float const dx(abs(x - xc)), dy(abs(y - yc)), dist(sqrt(dx*dx + dy*dy));
					if (dist > rval) continue;
					float const atten(0.6*dist*rval_inv);
					UNROLL_3X(data[y*stride + x].n[i_] *= atten;);
				}
			}
		}
	}

	void create_data(vector<vert_type_t> &data, vector<unsigned short> indices[NUM_LODS]) {
		//RESET_TIME;
		assert(zvals.size() == zvsize*zvsize);
		data.resize(stride*stride);
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
				data[y*stride + x] = vert_norm_comp(v, get_norm(ix)*light_scale);
			}
		}
		if (trees_enabled()) {
			init_tree_draw();
			apply_tree_ao_shadows(data);
		}
		for (unsigned i = 0; i < NUM_LODS; ++i) {
			indices[i].resize(4*(size>>i)*(size>>i));
			unsigned const step(1 << i);

			for (unsigned y = 0; y < size; y += step) {
				for (unsigned x = 0; x < size; x += step) {
					unsigned const vix(y*stride + x), iix(4*((y>>i)*(size>>i) + (x>>i)));
					indices[i][iix+0] = vix;
					indices[i][iix+1] = vix + step*stride;
					indices[i][iix+2] = vix + step*stride + step;
					indices[i][iix+3] = vix + step;
				}
			}
		}
		//PRINT_TIME("Create Data");
	}

	void create_texture() {
		assert(tid == 0);
		assert(!island);
		assert(zvals.size() == zvsize*zvsize);
		//RESET_TIME;
		unsigned char *data(new unsigned char[4*size*size]); // RGBA
		float const MESH_NOISE_SCALE = 0.003;
		float const MESH_NOISE_FREQ  = 80.0;
		float const dz_inv(1.0/(zmax - zmin)), noise_scale(MESH_NOISE_SCALE*mesh_scale_z);
		int k1, k2, k3, k4, dirt_tex_ix(-1), rock_tex_ix(-1);
		float t;

		for (unsigned i = 0; i < NTEX_DIRT; ++i) {
			if (lttex_dirt[i].id == DIRT_TEX) dirt_tex_ix = i;
			if (lttex_dirt[i].id == ROCK_TEX) rock_tex_ix = i;
		}
		assert(dirt_tex_ix >= 0 && rock_tex_ix >= 0);
		if (world_mode == WMODE_INF_TERRAIN) {create_xy_arrays(zvsize, MESH_NOISE_FREQ);}

		//#pragma omp parallel for schedule(static,1)
		for (unsigned y = 0; y < size; ++y) {
			for (unsigned x = 0; x < size; ++x) {
				float weights[NTEX_DIRT] = {0};
				unsigned const ix(y*zvsize + x);
				float const mh00(zvals[ix]), mh01(zvals[ix+1]), mh10(zvals[ix+zvsize]), mh11(zvals[ix+zvsize+1]);
				float const rand_offset((world_mode == WMODE_INF_TERRAIN) ? noise_scale*fast_eval_from_index(x, y, 0) : 0.0);
				float const relh1(relh_adj_tex + (min(min(mh00, mh01), min(mh10, mh11)) - zmin)*dz_inv + rand_offset);
				float const relh2(relh_adj_tex + (max(max(mh00, mh01), max(mh10, mh11)) - zmin)*dz_inv + rand_offset);
				get_tids(relh1, NTEX_DIRT-1, h_dirt, k1, k2, t);
				get_tids(relh2, NTEX_DIRT-1, h_dirt, k3, k4, t);
				bool const same_tid(k1 == k4);
				k2 = k4;
				
				if (same_tid) {
					t = 0.0;
				}
				else {
					float const relh(relh_adj_tex + (mh00 - zmin)*dz_inv);
					get_tids(relh, NTEX_DIRT-1, h_dirt, k1, k2, t);
				}
				float const sthresh[2] = {0.45, 0.7};
				float const vnz(get_norm(ix).z);
				float weight_scale(1.0);

				if (vnz < sthresh[1]) { // handle steep slopes (dirt/rock texture replaces grass texture)
					int const id(lttex_dirt[k1].id), id2(lttex_dirt[k2].id);

					if (id == GROUND_TEX || id2 == GROUND_TEX) { // ground/grass
						float const rock_weight((id == GROUND_TEX || id2 == ROCK_TEX) ? t : 0.0);
						weight_scale = CLIP_TO_01((vnz - sthresh[0])/(sthresh[1] - sthresh[0]));
						weights[rock_tex_ix] += (1.0 - weight_scale)*rock_weight;
						weights[dirt_tex_ix] += (1.0 - weight_scale)*(1.0 - rock_weight);
					}
					else if (id2 == SNOW_TEX) { // snow
						weight_scale = CLIP_TO_01(2.0f*(vnz - sthresh[0])/(sthresh[1] - sthresh[0]));
						weights[rock_tex_ix] += 1.0 - weight_scale;
					}
				}
				weights[k2] += weight_scale*t;
				weights[k1] += weight_scale*(1.0 - t);
				unsigned const off(4*(y*size + x));

				for (unsigned i = 0; i < NTEX_DIRT-1; ++i) { // Note: weights should sum to 1.0, so we can calculate w4 as 1.0-w0-w1-w2-w3
					data[off+i] = (unsigned char)(255.0*CLIP_TO_01(weights[i]));
				}
			} // for x
		} // for y
		setup_texture(tid, GL_MODULATE, 0, 0, 0, 0, 0);
		assert(tid > 0 && glIsTexture(tid));
		glTexImage2D(GL_TEXTURE_2D, 0, 4, size, size, 0, GL_RGBA, GL_UNSIGNED_BYTE, data); // internal_format = GL_COMPRESSED_RGBA - too slow
		glDisable(GL_TEXTURE_2D);
		delete [] data;
		//PRINT_TIME("Texture Upload");
	}

	float get_rel_dist_to_camera() const {
		return (p2p_dist_xy(get_camera_pos(), get_center()) - radius)/(get_tile_radius()*(X_SCENE_SIZE + Y_SCENE_SIZE));
	}
	float get_dist_to_camera_in_tiles() const {return get_rel_dist_to_camera()*get_tile_radius();}

	bool update_range() { // if returns 0, tile will be deleted
		float const dist(get_rel_dist_to_camera());
		if (dist > CLEAR_DIST_TILES) clear_vbo_tid(1,1);
		return (dist < DELETE_DIST_TILES);
	}

	bool is_visible() const {
		return camera_pdu.sphere_and_cube_visible_test(get_center(), radius, get_cube());
	}
	bool trees_are_distant() const {
		return (get_dist_to_camera_in_tiles() >= max(1.0, 5.0*calc_tree_size()));
	}
	bool use_low_tree_detail() const {return (ENABLE_TREE_LOD && trees_are_distant());}

	void init_tree_draw() {
		if (trees.generated) return; // already generate
		init_tree_dxoff = -xoff2;
		init_tree_dyoff = -yoff2;
		trees.gen_trees(x1+init_tree_dxoff, y1+init_tree_dyoff, x2+init_tree_dxoff, y2+init_tree_dyoff);
	}

	void update_tree_draw() {
		unsigned const desired_tlod(use_low_tree_detail() ? 1 : 2);

		if (tree_lod_level != desired_tlod) {
			tree_lod_level = desired_tlod;
			trees.clear_vbo_manager();
			trees.finalize(tree_lod_level == 1);
		}
		trees.vbo_manager.upload();
	}

	unsigned num_trees() const {return trees.size();}

	void draw_trees(vector<point> &trunk_pts, bool draw_branches, bool draw_leaves) const {
		glPushMatrix();
		bool const distant(trees_are_distant());
		vector3d const xlate(((xoff - xoff2) - init_tree_dxoff)*DX_VAL, ((yoff - yoff2) - init_tree_dyoff)*DY_VAL, 0.0);
		translate_to(xlate);
		
		if (draw_branches) {
			if (distant) {
				trees.add_trunk_pts(xlate, trunk_pts);
			}
			else {
				trees.draw_branches(0, xlate, &trunk_pts);
			}
		}
		if (draw_leaves) {
			bool const cull(ENABLE_TREE_BFC && distant);
			bool const draw_all(use_low_tree_detail() || camera_pdu.sphere_completely_visible_test(get_center(), radius));
			if (cull) {glEnable (GL_CULL_FACE);}
			trees.draw_leaves(0, draw_all, xlate);
			if (cull) {glDisable(GL_CULL_FACE);}
		}
		glPopMatrix();
	}

	void draw(vector<vert_type_t> &data, vector<unsigned short> indices[NUM_LODS], bool reflection_pass) {
		assert(size > 0);

		if (!DISABLE_TEXTURES) {
			if (tid == 0) create_texture();
			
			if (tid > 0) {
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, tid);
			}
		}
		glPushMatrix();
		glTranslatef(((xoff - xoff2) - init_dxoff)*DX_VAL, ((yoff - yoff2) - init_dyoff)*DY_VAL, 0.0);
		if (tid > 0) set_landscape_texgen(1.0, (-x1 - init_dxoff), (-y1 - init_dyoff), MESH_X_SIZE, MESH_Y_SIZE);
		unsigned const ptr_stride(sizeof(vert_type_t));

		if (vbo == 0) {
			create_data(data, indices);
			vbo = create_vbo();
			bind_vbo(vbo, 0);
			upload_vbo_data(&data.front(), data.size()*ptr_stride, 0);

			for (unsigned i = 0; i < NUM_LODS; ++i) {
				assert(ivbo[i] == 0);
				ivbo[i] = create_vbo();
				bind_vbo(ivbo[i], 1);
				assert(!indices[i].empty());
				upload_vbo_data(&(indices[i].front()), indices[i].size()*sizeof(unsigned short), 1);
			}
		}
		unsigned lod_level(reflection_pass ? min(NUM_LODS-1, 1U) : 0);
		float dist(get_dist_to_camera_in_tiles());

        while (dist > (reflection_pass ? 1.0 : 2.0) && lod_level+1 < NUM_LODS) {
            dist /= 2;
            ++lod_level;
        }
		assert(lod_level < NUM_LODS);
		assert(vbo > 0 && ivbo[lod_level] > 0);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo[lod_level], 1);
		unsigned const isz(size >> lod_level);
		
		// can store normals in a normal map texture, but a vertex texture fetch is slow
		glVertexPointer(3, GL_FLOAT, ptr_stride, 0);
		glNormalPointer(GL_BYTE, ptr_stride, (void *)sizeof(point)); // was GL_FLOAT
		glDrawRangeElements(GL_QUADS, 0, (unsigned)data.size(), 4*isz*isz, GL_UNSIGNED_SHORT, 0); // requires GL/glew.h
		bind_vbo(0, 0);
		bind_vbo(0, 1);
		glPopMatrix();
		if (tid > 0) disable_textures_texgen();
	}
}; // tile_t


class tile_draw_t {

	typedef map<tile_xy_pair, tile_t*> tile_map;
	tile_map tiles;
	vector<point> tree_trunk_pts;

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
		tree_trunk_pts.clear();
	}

	void update() {
		//RESET_TIME;
		assert(MESH_X_SIZE == MESH_Y_SIZE); // limitation, for now
		point const camera(get_camera_pos() - point((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
		int const tile_radius(int(1.5*get_tile_radius()) + 1);
		int const toffx(int(0.5*camera.x/X_SCENE_SIZE)), toffy(int(0.5*camera.y/Y_SCENE_SIZE));
		int const x1(-tile_radius + toffx), y1(-tile_radius + toffy);
		int const x2( tile_radius + toffx), y2( tile_radius + toffy);
		unsigned const init_tiles((unsigned)tiles.size());
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
				if (tile.get_rel_dist_to_camera() >= CREATE_DIST_TILES) continue; // too far away to create
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


	static void setup_terrain_textures(shader_t &s, unsigned start_tu_id, bool use_sand) {
		unsigned const base_tsize(get_norm_texels());

		for (int i = 0; i < (use_sand ? NTEX_SAND : NTEX_DIRT); ++i) {
			int const tid(use_sand ? lttex_sand[i].id : lttex_dirt[i].id);
			float const tscale(float(base_tsize)/float(get_texture_size(tid, 0))); // assumes textures are square
			float const cscale((world_mode == WMODE_INF_TERRAIN && tid == GROUND_TEX) ? 0.5 : 1.0); // darker grass
			unsigned const tu_id(start_tu_id + i);
			select_multitex(tid, tu_id, 0);
			std::ostringstream oss1, oss2, oss3;
			oss1 << "tex" << tu_id;
			oss2 << "ts"  << tu_id;
			oss3 << "cs"  << tu_id;
			s.add_uniform_int(  oss1.str().c_str(), tu_id);
			s.add_uniform_float(oss2.str().c_str(), tscale);
			s.add_uniform_float(oss3.str().c_str(), cscale);
		}
		set_multitex(0);
	}

	static void setup_mesh_draw_shaders(shader_t &s, float wpz, bool blend_textures_in_shaders, bool reflection_pass) {
		s.setup_enabled_lights();
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		s.set_vert_shader("texture_gen.part+tiled_mesh");
		s.set_frag_shader(blend_textures_in_shaders ? "linear_fog.part+tiled_mesh" : "linear_fog.part+multitex_2");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		s.add_uniform_int("tex1", 1);
		s.add_uniform_float("water_plane_z", (is_water_enabled() ? wpz : zmin));
		s.add_uniform_float("water_atten", WATER_COL_ATTEN*mesh_scale);
		s.add_uniform_float("normal_z_scale", (reflection_pass ? -1.0 : 1.0));
		if (blend_textures_in_shaders) setup_terrain_textures(s, 2, 0);
	}

	float draw(bool add_hole, float wpz, bool reflection_pass) {
		float zmin(FAR_CLIP);
		glDisable(GL_NORMALIZE);
		set_array_client_state(1, 0, 1, 0);
		unsigned num_drawn(0), num_trees(0);
		unsigned long long mem(0);
		vector<tile_t::vert_type_t> data;
		vector<unsigned short> indices[NUM_LODS];
		static point last_sun(all_zeros), last_moon(all_zeros);

		if (sun_pos != last_sun || moon_pos != last_moon) {
			clear_vbos_tids(1,0); // light source changed, clear vbos and build new shadow map
			last_sun  = sun_pos;
			last_moon = moon_pos;
		}
		shader_t s;
		setup_mesh_draw_shaders(s, wpz, 1, reflection_pass);
		
		if (world_mode == WMODE_INF_TERRAIN && show_fog && is_water_enabled() && !reflection_pass) {
			draw_water_edge(wpz); // Note: doesn't take into account waves
		}
		setup_mesh_lighting();
		vector<pair<float, tile_t *> > to_draw;

		for (tile_map::const_iterator i = tiles.begin(); i != tiles.end(); ++i) {
			assert(i->second);
			if (DEBUG_TILES) mem += i->second->get_gpu_memory();
			if (add_hole && i->first.x == 0 && i->first.y == 0) continue;
			float const dist(i->second->get_rel_dist_to_camera());
			if (dist > DRAW_DIST_TILES) continue; // too far to draw
			zmin = min(zmin, i->second->get_zmin());
			if (!i->second->is_visible()) continue;
			to_draw.push_back(make_pair(dist, i->second));
		}
		sort(to_draw.begin(), to_draw.end()); // sort front to back to improve draw time through depth culling

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			num_trees += to_draw[i].second->num_trees();
			to_draw[i].second->draw(data, indices, reflection_pass);
		}
		s.end_shader();
		if (DEBUG_TILES) cout << "tiles drawn: " << to_draw.size() << " of " << tiles.size() << ", trees: " << num_trees << ", gpu mem: " << mem/1024/1024 << endl;
		run_post_mesh_draw();
		if (trees_enabled()) {draw_trees(to_draw, reflection_pass);}
		return zmin;
	}

	void draw_trees(vector<pair<float, tile_t *> > const &to_draw, bool reflection_pass) {
		for (unsigned i = 0; i < to_draw.size(); ++i) {
			if (trees_enabled()) {to_draw[i].second->update_tree_draw();}
		}
		shader_t s;
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		colorRGBA const orig_fog_color(setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0));
		s.add_uniform_float("tex_scale_t", 5.0);
		set_color(get_tree_trunk_color(T_PINE, 0)); // all a constant color

		for (unsigned i = 0; i < to_draw.size(); ++i) { // branches
			to_draw[i].second->draw_trees(tree_trunk_pts, 1, 0);
		}
		s.add_uniform_float("tex_scale_t", 1.0);
		s.end_shader();
		s.set_prefix("#define USE_LIGHT_COLORS",  0); // VS
		s.set_prefix("#define USE_GOOD_SPECULAR", 1); // VS
		s.set_prefix("#define USE_QUADRATIC_FOG", 1); // FS
		s.setup_enabled_lights(2);
		s.set_vert_shader("ads_lighting.part*+pine_tree");
		s.set_frag_shader("linear_fog.part+pine_tree");
		//s.set_geom_shader("pine_tree", GL_POINTS, GL_TRIANGLE_STRIP, 120); // actually outputs quads
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		s.add_uniform_float("min_alpha", 0.75);
		check_gl_error(302);
		
		if (!tree_trunk_pts.empty()) { // color/texture already set above
			assert(!(tree_trunk_pts.size() & 1));
			s.add_uniform_float("camera_facing_scale", 1.0);
			select_texture(WHITE_TEX, 0); // enable=0
			get_tree_trunk_color(T_PINE, 1).do_glColor();
			zero_vector.do_glNormal();
			set_array_client_state(1, 0, 0, 0);
			glVertexPointer(3, GL_FLOAT, sizeof(point), &tree_trunk_pts.front());
			glDrawArrays(GL_LINES, 0, (unsigned)tree_trunk_pts.size());
			tree_trunk_pts.resize(0);
		}
		set_specular(0.2, 8.0);

		for (unsigned i = 0; i < to_draw.size(); ++i) { // leaves
			s.add_uniform_float("camera_facing_scale", (to_draw[i].second->use_low_tree_detail() ? 1.0 : 0.0));
			to_draw[i].second->draw_trees(tree_trunk_pts, 0, 1);
		}
		assert(tree_trunk_pts.empty());
		set_specular(0.0, 1.0);
		s.end_shader();
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
}; // tile_draw_t


tile_draw_t terrain_tile_draw;


tile_t *get_tile_from_xy(tile_xy_pair const &tp) {

	return terrain_tile_draw.get_tile_from_xy(tp);
}


void draw_vert_color(colorRGBA c, float x, float y, float z) {

	c.do_glColor();
	glVertex3f(x, y, z);
}


void fill_gap(float wpz, bool reflection_pass) {

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
	shader_t s;
	terrain_tile_draw.setup_mesh_draw_shaders(s, wpz, 0, reflection_pass);

	// draw +x
	build_xy_mesh_arrays(&xv.front(), &yv[MESH_Y_SIZE], MESH_X_SIZE, 1);
	glBegin(GL_QUAD_STRIP);

	for (int x = 0; x < MESH_X_SIZE; ++x) {
		vertex_normals[MESH_Y_SIZE-1][x].do_glNormal();
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE-DY_VAL, mesh_height[MESH_Y_SIZE-1][x]);
		draw_vert_color(color, get_xval(x), Y_SCENE_SIZE       , fast_eval_from_index(x, 0));
	}
	glEnd();
	float const z_end_val(fast_eval_from_index(MESH_X_SIZE-1, 0));

	// draw +y
	build_xy_mesh_arrays(&xv[MESH_X_SIZE], &yv.front(), 1, MESH_Y_SIZE+1);
	glBegin(GL_QUAD_STRIP);

	for (int y = 0; y < MESH_X_SIZE; ++y) {
		vertex_normals[y][MESH_X_SIZE-1].do_glNormal();
		draw_vert_color(color, X_SCENE_SIZE-DX_VAL, get_yval(y), mesh_height[y][MESH_X_SIZE-1]);
		draw_vert_color(color, X_SCENE_SIZE       , get_yval(y), fast_eval_from_index(0, y));
	}

	// draw corner quad
	draw_vert_color(color, X_SCENE_SIZE-DX_VAL, Y_SCENE_SIZE, z_end_val);
	draw_vert_color(color, X_SCENE_SIZE       , Y_SCENE_SIZE, fast_eval_from_index(0, MESH_Y_SIZE));
	glEnd();

	s.end_shader();
	disable_textures_texgen();
	glDisable(GL_COLOR_MATERIAL);
	run_post_mesh_draw();
	//PRINT_TIME("Fill Gap");
}


float draw_tiled_terrain(bool add_hole, float wpz, bool reflection_pass) {

	//RESET_TIME;
	bool const vbo_supported(setup_gen_buffers());
		
	if (!vbo_supported) {
		cout << "Warning: VBOs not supported, so tiled mesh cannot be drawn." << endl;
		return zmin;
	}
	terrain_tile_draw.update();
	float const zmin(terrain_tile_draw.draw(add_hole, wpz, reflection_pass));
	if (add_hole) fill_gap(wpz, reflection_pass); // need to fill the gap on +x/+y
	//glFinish(); PRINT_TIME("Tiled Terrain Draw"); //exit(0);
	return zmin;
}


void clear_tiled_terrain() {
	terrain_tile_draw.clear();
}

void reset_tiled_terrain_state() {
	terrain_tile_draw.clear_vbos_tids(1,1);
}


