// 3D World - Snow Accumulation and Renderign code
// by Frank Gennari
// 5/5/10
#include "3DWorld.h"
#include "mesh.h"
#include "gl_ext_arb.h"
#include "shaders.h"


unsigned const VOXELS_PER_DIV = 8; // 1024 for 128 vertex mesh
unsigned const MAX_STRIP_LEN  = 200; // larger = faster, less overhead; smaller = smaller edge strips, better culling
int      const Z_CHECK_RANGE  = 1; // larger = smoother and fewer strips, but longer preprocessing
bool     const USE_VBOS       = 1; // faster drawing but more GPU resources

bool has_snow(0);
point vox_delta;
map<int, unsigned> x_strip_map;

extern bool disable_shaders;
extern int display_mode, read_snow_file, write_snow_file;
extern unsigned num_snowflakes;
extern float ztop, zbottom, temperature, snow_depth, snow_random;
extern vector3d cview_dir;
extern char *snow_file;
extern vector<coll_obj> coll_objects;


typedef short coord_type;
typedef unsigned short count_type;
unsigned const MAX_COUNT((1 << (sizeof(count_type) << 3)) - 1);

void increment_printed_number(unsigned num);


// **************** LOW-LEVEL CONTAINER CLASSES ****************

struct zval_avg {
	count_type c;
	float z;
	zval_avg(count_type c_=0, float z_=0.0) : c(c_), z(z_) {}
	// Note: With a 1024x1024 voxel grid we should get an average of 1 count per 1M snow
	//       so we can have at max 64M snowflakes.
	//       However, we can get snow to stack up at a vertical edge so we need to clamp the count
	void update(float zval) {if (c < MAX_COUNT) {++c; z += zval;}}
	bool valid() const {return (c > 0);}
	float getz() const {return z/c;}
};


struct voxel_t {
	coord_type p[3];

	voxel_t(void) {}
	voxel_t(coord_type x, coord_type y, coord_type z) {p[0] = x; p[1] = y; p[2] = z;}
	
	voxel_t(point const &pt) {
		p[0] = int(vox_delta.x*(pt.x + X_SCENE_SIZE) + 0.5);
		p[1] = int(vox_delta.y*(pt.y + Y_SCENE_SIZE) + 0.5);
		p[2] = int(vox_delta.z*(pt.z - czmin));
	}
	point get_pt() const {
		return point((-X_SCENE_SIZE + p[0]/vox_delta.x), (-Y_SCENE_SIZE + p[1]/vox_delta.y), (czmin + p[2]/vox_delta.z));
	}
	bool operator<(voxel_t const &v) const { // sorted {x, y, z}
		if (p[0] < v.p[0]) return 1;
		if (p[0] > v.p[0]) return 0;
		if (p[1] < v.p[1]) return 1;
		if (p[1] > v.p[1]) return 0;
		return (p[2] < v.p[2]);
	}
	// unused, but useful if implemented with a hash_map
	size_t hash() const {return (345379*size_t(p[0]) + 988439*size_t(p[1]) + 632447*size_t(p[2]));}
	bool operator==(voxel_t const &v) const {return (p[0] == v.p[0] && p[1] == v.p[1] && p[2] == v.p[2]);}
};


struct voxel_z_pair {
	voxel_t v;
	zval_avg z;
	voxel_z_pair() {}
	voxel_z_pair(voxel_t const &v_, zval_avg const &z_=zval_avg()) : v(v_), z(z_) {}
};


struct strip_entry {
	point p;
	vector3d n;
	strip_entry() {}
	strip_entry(point const &p_) : p(p_) {}
};

class strip_vect_t : public vector<strip_entry> {
	float zval;
	bool z_valid;

public:
	strip_vect_t() : zval(0.0), z_valid(0) {}

	void add(voxel_z_pair const &vz, float d) {
		point p(vz.v.get_pt());
		
		if (vz.z.valid()) {
			zval  = vz.z.getz();
			p.z   = zval + vz.z.c*d; // zval is average zval + snow depth
			zval -= 0.1/vox_delta.z; // shift to move sligltly under the object (optional)

			if (!z_valid) {
				for (unsigned i = 0; i < size(); ++i) {
					operator[](i).p.z = zval;
				}
				z_valid = 1;
			}
		}
		else if (z_valid) {
			p.z = zval;
		}
		push_back(strip_entry(p));
	}
	void reset() {
		resize(0);
		z_valid = 0;
	}
};


// **************** BLOCK ALLOCATION ****************

unsigned const BLOCK_SZ = 65536; // 1.5MB blocks

class strip_block_alloc {

	vector<strip_entry *> blocks;
	unsigned pos;

public:
	strip_block_alloc() : pos(0) {}
	unsigned get_pos() const {return pos;}

	~strip_block_alloc() {
		for (vector<strip_entry *>::iterator i = blocks.begin(); i != blocks.end(); ++i) {
			delete [] *i;
		}
	}

	strip_entry *alloc(unsigned sz) {
		assert(sz <= 2*BLOCK_SZ); // could use larger block size or singles vector

		if (blocks.empty() || (pos + sz) > BLOCK_SZ) {
			blocks.push_back(new strip_entry[BLOCK_SZ]);
			pos = 0;
		}
		strip_entry *const cur(blocks.back() + pos);
		pos += sz;
		return cur;
	}

	void status() const {
		cout << "blocks: " << blocks.size() << ", pos: " << pos
			 << ", mem: "  << blocks.size()*BLOCK_SZ*sizeof(strip_entry) << endl;
	}

	unsigned get_cur_block_id() const {
		assert(!blocks.empty());
		return (blocks.size() - 1);
	}
};

strip_block_alloc sb_alloc;


// **************** CORE CLASSES ****************

class snow_renderer; // forward declaration

class strip_t {

	strip_entry *strips; // block allocated
	unsigned size, block_id, block_pos;

public:
	friend class snow_renderer;
	bool is_edge;
	int xval, y_start;

	strip_t(bool edge=0) : strips(NULL), size(0), block_id(0), block_pos(0), is_edge(edge), xval(0), y_start(0) {}
	unsigned get_size() const {return size;}

	point get_pos(unsigned pos) const {
		assert(strips && pos < size);
		return strips[pos].p;
	}
	point get_norm(unsigned pos) const {
		assert(strips && pos < size);
		return strips[pos].n;
	}

	void finalize(strip_vect_t const &strip_vect, int xval_, int y_start_) { // can only call this once
		xval    = xval_;
		y_start = y_start_;

		// allocate data
		assert(size == 0);
		size      = strip_vect.size();
		strips    = sb_alloc.alloc(size);
		block_id  = sb_alloc.get_cur_block_id();
		block_pos = sb_alloc.get_pos();
		assert(strips);
		assert(size <= block_pos);
		block_pos -= size; // get starting position

		for (unsigned i = 0; i < size; ++i) {
			strips[i] = strip_vect[i];
		}

		// calculate normals
		assert(size >= 4);
		vector<vector3d> nt(size-2); // triangle face normals

		for (unsigned i = 0; i < nt.size(); ++i) { // calculate triangle normals
			nt[i] = cross_product((strips[i].p - strips[i+1].p), (strips[i+2].p - strips[i+1].p)).get_norm();
			if (nt[i].z < 0.0) nt[i].negate(); // normal is always +z
		}
		for (unsigned i = 0; i < size; ++i) { // average triangle normals to get vertex normals
			unsigned const k1((i > 2) ? i-2 : 0), k2(min(nt.size()-1, i));
			assert(k2 >= k1);
			vector3d &n(strips[i].n);
			n = zero_vector;
			for (unsigned k = k1; k <= k2; ++k) {n += nt[k];}
			n /= (k2 - k1 + 1); // average triangle normals
		}
	}
};


class voxel_map : public map<voxel_t, zval_avg> { // hash_map or Google sparse_hashmap?
public:
	zval_avg find_adj_z(voxel_t &v, zval_avg const &zv_old, float depth, voxel_map *cur_x_map=NULL);
	bool read(char const *const fn);
	bool write(char const *const fn) const;
};


struct data_block { // packed voxel_z_pair for read/write
	coord_type p[3];
	count_type c;
	float z;

	void add_to_map(voxel_map &vmap) const {
		voxel_t v;
		for (unsigned i = 0; i < 3; ++i) v.p[i] = p[i];
		vmap[v] = zval_avg(c, z);
	}
	void set_from_map_iter(voxel_map::const_iterator it) {
		for (unsigned i = 0; i < 3; ++i) p[i] = it->first.p[i];
		c = it->second.c;
		z = it->second.z;
	}
};


// this tends to take a large fraction of the preprocessing time
zval_avg voxel_map::find_adj_z(voxel_t &v, zval_avg const &zv_old, float depth, voxel_map *cur_x_map) {

	coord_type best_dz(0);
	zval_avg res;

	for (coord_type dz = -Z_CHECK_RANGE; dz <= Z_CHECK_RANGE; ++dz) { // check adjacent z values (at, above, below)
		voxel_t v2(v);
		v2.p[2] += dz;
		const_iterator it(find(v2));
		if (it == end()) continue; // zero entry
		zval_avg const z2(it->second);
		assert(z2.valid());
		if (zv_old.valid() && fabs(z2.getz() - zv_old.getz()) > depth) continue; // delta z too large
		
		if (cur_x_map) {
			erase(v2);
			(*cur_x_map)[v2] = z2;
		}
		if (!res.valid() || abs(dz) < abs(best_dz)) best_dz = dz;
		res.c += z2.c;
		res.z += z2.z;
	}
	if (res.valid()) v.p[2] += best_dz;
	return res;
}


bool voxel_map::read(char const *const fn) {

	FILE *fp;
	assert(fn != NULL);
	if (!open_file(fp, fn, "snow map", "rb")) return 0;
	cout << "Reading snow file from " << fn << endl;
	size_t const n(fread(&vox_delta, sizeof(float), 3, fp));
	assert(n == 3);
	unsigned map_size(0);
	size_t const sz_read(fread(&map_size, sizeof(unsigned), 1, fp));
	assert(sz_read == 1);
	
	for (unsigned i = 0; i < map_size; ++i) {
		data_block data;
		size_t const n(fread(&data, sizeof(data_block), 1, fp));
		assert(n == 1);
		data.add_to_map(*this);
	}
	fclose(fp);
	return 1;
}


bool voxel_map::write(char const *const fn) const {

	FILE *fp;
	assert(fn != NULL);
	if (!open_file(fp, fn, "snow map", "wb")) return 0;
	cout << "Writing lighting file to " << fn << endl;
	size_t const n(fwrite(&vox_delta, sizeof(float), 3, fp));
	assert(n == 3);
	unsigned const map_size(size());
	size_t const sz_write(fwrite(&map_size, sizeof(unsigned), 1, fp));
	assert(sz_write == 1);

	for (const_iterator i = begin(); i != end(); ++i) {
		data_block data;
		data.set_from_map_iter(i);
		size_t const n(fwrite(&data, sizeof(data_block), 1, fp));
		assert(n == 1);
	}
	fclose(fp);
	return 1;
}


// **************** VBO/ELEMENT INDEX CODE ****************


class snow_renderer {

	unsigned vbo, ivbo;
	float last_x;
	vector<vert_norm> data;
	vector<unsigned> indices, strip_offsets;
	map<point, unsigned> vmap[2]; // {prev, next} rows

public:
	snow_renderer() : vbo(0), ivbo(0), last_x(0.0) {}
	~snow_renderer() {free_vbos();}
	bool empty() const {return data.empty();}

	void add_all_strips(vector<strip_t> const &strips) {
		unsigned nquads(0);

		for (vector<strip_t>::const_iterator i = strips.begin(); i != strips.end(); ++i) {
			nquads += (i->get_size() - 2)/2;
		}
		indices.reserve(4*nquads);
		data.reserve(5*nquads/4); // 20% extra
		strip_offsets.reserve(strips.size()+1);

		for (vector<strip_t>::const_iterator i = strips.begin(); i != strips.end(); ++i) {
			strip_offsets.push_back(indices.size());
			add_strip(*i);
		}
		strip_offsets.push_back(indices.size());
		assert(indices.size() == 4*nquads);
	}

	void add_strip(strip_t const &s) {
		unsigned const size(s.get_size());
		assert(size >= 4 && !(size & 1)); // must be even
		float const xval(s.strips[0].p.x);

		if (xval != last_x) {
			vmap[0].clear();
			vmap[0].swap(vmap[1]);
			last_x = xval;
		}
		for (unsigned i = 0; i+2 < size; i += 2) { // iterate as a quad strip
			add(s.strips[i+0].p, s.strips[i+0].n, 0);
			add(s.strips[i+1].p, s.strips[i+1].n, 1);
			add(s.strips[i+3].p, s.strips[i+3].n, 1);
			add(s.strips[i+2].p, s.strips[i+2].n, 0);
		}
	}

private:
	unsigned add(point const &v, vector3d const &n, unsigned map_ix) { // can't be called after finalize()
		assert(vbo == 0 && ivbo == 0);
		map<point, unsigned>::const_iterator it(vmap[map_ix].find(v));
		unsigned ix(0);

		if (it != vmap[map_ix].end()) { // existing vertex
			ix = it->second;
			assert(ix < data.size());
			data[ix].n = (data[ix].n + n)*0.5; // average the normals???
		}
		else {
			ix = data.size();
			vmap[map_ix][v] = ix;
			data.push_back(vert_norm(v, n));
		}
		indices.push_back(ix);
		return ix;
	}

	void calc_shadows() {
		int index(-1);
		point const lpos(get_light_pos());

		for (vector<vert_norm>::iterator i = data.begin(); i != data.end(); ++i) {
			i->n.normalize();
			bool shadowed(0);

			if (index >= 0 && coll_objects[index].line_intersect(i->v, lpos)) {
				shadowed = 1;
			}
			else {
				index = -1;
				shadowed = !is_visible_to_light_cobj(i->v, get_light(), 0.0, -1, 1, &index);
				if (index >= 0) assert(index < (int)coll_objects.size());
			}
			if (shadowed) i->n *= 0.001; // scale to a small (zero-ish) value
		}
	}

	void upload_vbo() {
		assert(!data.empty());
		if (vbo == 0) vbo = create_vbo();
		bind_vbo(vbo, 0);
		upload_vbo_data(&data.front(), data.size()*sizeof(vert_norm), 0);
		bind_vbo(0, 0);
	}

	void upload_ivbo() {
		assert(!indices.empty());
		if (ivbo == 0) ivbo = create_vbo();
		bind_vbo(ivbo, 1);
		upload_vbo_data(&indices.front(), indices.size()*sizeof(unsigned), 1);
		bind_vbo(0, 1);
	}

public:
	void free_vbos() {
		delete_vbo(vbo);
		delete_vbo(ivbo);
		vbo = ivbo = 0;
	}

	void update_shadows() {
		calc_shadows();
		upload_vbo();
	}

	void update_region(unsigned strip_ix, unsigned strip_pos, unsigned strip_len, float new_z) {
		// FIXME: update a range at a time?
		assert(vbo);
		bind_vbo(vbo, 0);
		assert(strip_ix+1 < strip_offsets.size());
		assert(strip_len >= 4); // at least one quad
		unsigned const cur_six(strip_offsets[strip_ix]), next_six(strip_offsets[strip_ix+1]);
		unsigned const num_quads((strip_len-2)/2), quad_ix(min(strip_pos/2, num_quads-1));
		assert((next_six - cur_six) == 4*num_quads); // error check
		unsigned const start_index_ix(cur_six + 4*quad_ix); // quad vertex index

		for (unsigned i = 0; i < 4; ++i) { // 4 points on the quad
			unsigned index_ix(start_index_ix + i);
			assert(index_ix < indices.size());
			assert(index_ix < next_six);
			unsigned const data_ix(indices[index_ix]);
			assert(data_ix < data.size());
			vector3d const norm(((i < 2) ? 1.0 : -1.0), ((i & 1) ? -1.0 : 1.0), 0.0);
			data[data_ix].n   = (data[data_ix].n + norm.get_norm())*0.5; // FIXME: very approximate
			data[data_ix].v.z = new_z;
			upload_vbo_sub_data(&data[data_ix], data_ix*sizeof(vert_norm), sizeof(vert_norm), 0);
		}
		bind_vbo(0, 0);
	}

	void finalize() { // can only be called once
		assert(vbo == 0 && ivbo == 0);
		assert((indices.size() & 3) == 0); // must be a multiple of 4
		upload_ivbo();
		vmap[0].clear();
		vmap[1].clear();
	}

	void draw() {
		if (vbo  == 0) upload_vbo ();
		if (ivbo == 0) upload_ivbo();
		assert(vbo != 0 && ivbo != 0);
		bind_vbo(vbo,  0);
		bind_vbo(ivbo, 1);
		glVertexPointer(3, GL_FLOAT, sizeof(vert_norm), 0);
		glNormalPointer(   GL_FLOAT, sizeof(vert_norm), (void *)sizeof(point));
		glDrawRangeElements(GL_QUADS, 0, data.size(), indices.size(), GL_UNSIGNED_INT, 0); // requires GL/glew.h
		//glDrawElements(GL_QUADS, indices.size(), GL_UNSIGNED_INT, 0);
		bind_vbo(0, 0);
		bind_vbo(0, 1);
	}

	void show_stats() const {
		unsigned const dmem(data.size()*sizeof(vert_norm));
		unsigned const imem(indices.size()*sizeof(unsigned));
		cout << "verts: " << data.size() << ", quads: " << indices.size()/4 << endl;
		cout << "mem: " << (dmem + imem) << ", vmem: " << (dmem*(vbo != 0) + imem*(ivbo != 0)) << endl;
	}
};

snow_renderer snow_draw;


// **************** TOP LEVEL FUNCTIONS ****************

vector<strip_t> snow_strips;


bool get_mesh_ice_pt(point const &p1, point &p2) {

	int const xpos(get_xpos(p1.x)), ypos(get_ypos(p1.y));
	if (point_outside_mesh(xpos, ypos)) return 0; // shouldn't get here (much)
	float const mh(interpolate_mesh_zval(p1.x, p1.y, 0.0, 0, 1));
	float const h(max(mh, water_matrix[ypos][xpos])); // mesh or ice (water)
	p2.assign(p1.x, p1.y, h);
	return 1;
}


vector3d get_rand_snow_vect(float amount) {

	return vector3d(amount*snow_random*rgauss(), amount*snow_random*rgauss(), 0.0);
}


void create_snow_map(voxel_map &vmap) {

	// distribute snowflakes over the scene and build the voxel map of hits
	int const num_per_dim(1024*(unsigned)sqrt((float)num_snowflakes)); // in M, so sqrt div by 1024
	float const zval(max(ztop, czmax));
	unsigned progress(0);
	cout << "Snow accumulation progress (out of " << num_per_dim << "): 0";

	#pragma omp parallel for schedule(static,1)
	for (int y = 0; y < num_per_dim; ++y) {
		#pragma omp critical(snow_prog_update)
		increment_printed_number(progress++);

		for (int x = 0; x < num_per_dim; ++x) {
			point pos1(-X_SCENE_SIZE + 2.0*X_SCENE_SIZE*x/num_per_dim,
				       -Y_SCENE_SIZE + 2.0*Y_SCENE_SIZE*y/num_per_dim, zval);

			// add slightly more randomness for numerical precision reasons
			for (unsigned d = 0; d < 2; ++d) {
				pos1[d] += SMALL_NUMBER*signed_rand_float();
			}
			point pos2;
			if (!get_mesh_ice_pt(pos1, pos2)) continue; // invalid point
			assert(pos2.z < pos1.z);
			pos1 += get_rand_snow_vect(1.0); // add some gaussian randomness for better distribution
			point cpos;
			vector3d cnorm;
			int cindex;
			
			while (check_coll_line_exact(pos1, pos2, cpos, cnorm, cindex, 0.0, -1, 0, 0, 1, 1)) {
				if (cnorm.z > 0.0) { // collision with a surface that points up
					pos2 = cpos;
					break;
				}
				// collision with vertical or bottom surface
				pos1 = cpos;
				
				if (!get_mesh_ice_pt(pos1, pos2)) { // invalid point
					pos2 = cpos;
					break;
				}
				float const val(CLIP_TO_01((pos1.z - zbottom)/(zval - zbottom)));
				vector3d delta(get_rand_snow_vect(val));
				if (dot_product(delta, cnorm) < 0.0) delta.negate();
				pos1   += delta; // push a random amount away from the object
				pos1.z -= SMALL_NUMBER; // ensure progress is made
			}
			#pragma omp critical(snow_map_update)
			vmap[voxel_t(pos2)].update(pos2.z);
		}
	}
	cout << endl;
}


void add_strip(strip_vect_t const &strip, bool is_edge, int xval, int y_start) {

	assert(!strip.empty());
	snow_strips.push_back(strip_t(is_edge));
	snow_strips.back().finalize(strip, xval, y_start);
}


void create_snow_strips(voxel_map &vmap) {

	// create strips of snow for rendering
	voxel_map cur_x_map, last_x_map;
	float const delta_depth(snow_depth*VOXELS_PER_DIV*VOXELS_PER_DIV*XY_MULT_SIZE/(1024.0*1024.0*num_snowflakes));
	unsigned n_strips(0), n_edge_strips(0), strip_len(0), edge_strip_len(0);
	int last_x(0);
	strip_vect_t strip, edge_strip;
	vector<voxel_z_pair> vs;

	while (!vmap.empty()) {
		pair<voxel_t, zval_avg> start(*vmap.begin());
		voxel_t v1(start.first);
		zval_avg zv(start.second);
		assert(zv.valid());

		if (v1.p[0] != last_x) { // we moved on to the next x-value, so update the x maps
			last_x = v1.p[0];
			last_x_map.clear();
			cur_x_map.swap(last_x_map);
			assert(x_strip_map.find(last_x) == x_strip_map.end()); // map should guarantee strictly increasing x
			x_strip_map[last_x] = snow_strips.size();
		}
		vmap.erase(v1);
		cur_x_map[v1] = zv;
		vs.resize(0);
		--v1.p[1];
		vs.push_back(voxel_z_pair(v1)); // zero start
		++v1.p[1];
		vs.push_back(voxel_z_pair(v1, zv));
		
		while (1) { // generate a strip in y with constant x
			++v1.p[1];
			zv = vmap.find_adj_z(v1, zv, snow_depth, &cur_x_map);
			//if (!zv.valid()) --v1.p[1]; // move back one step
			vs.push_back(voxel_z_pair(v1, zv));
			if (!zv.valid()) break; // end of strip
		}
		unsigned const sz(vs.size()), num_parts((sz - 1)/MAX_STRIP_LEN + 1); // ceiling

		for (unsigned n = 0; n < num_parts; ++n) {
			unsigned const start_pos(n*MAX_STRIP_LEN);
			unsigned end_pos(min(sz, start_pos+MAX_STRIP_LEN+1));
			unsigned num_ends(0), last_edge(0);
			int const y_start(vs[start_pos].v.p[1]);
			int edge_y_start(y_start);
			strip.reset();
			edge_strip.reset();

			if ((sz - end_pos) < MAX_STRIP_LEN/10) { // don't create tiny strips
				end_pos = sz;
				++n;
			}
			for (unsigned i = start_pos; i < end_pos; ++i) {
				bool const end_element(i == 0 || i+1 == sz);
				voxel_t v2(vs[i].v);
				++v2.p[0]; // move to next x row
				zval_avg z2(vmap.find_adj_z(v2, vs[i].z, snow_depth));
				if (end_element) z2.c = 0; // zero terminate start/end points
				strip.add(vs[i], delta_depth); // first edge
				strip.add(voxel_z_pair(v2, z2), delta_depth); // second edge

				// generate edge strips
				if ((end_pos - start_pos) <= 3) continue; // too small for edge srtips
				voxel_t v3(vs[i].v);
				--v3.p[0]; // move to prev x row
				zval_avg z3(last_x_map.find_adj_z(v3, vs[i].z, snow_depth));
				
				if (!end_element && !z3.valid()) {
					last_edge = edge_strip.size() + 2;
					++num_ends;
				}
				if (num_ends == 0) {
					edge_strip.reset();
					edge_y_start = v3.p[1];
				}
				if (end_element) z3.c = 0; // zero terminate start/end points
				edge_strip.add(vs[i], delta_depth); // first edge
				edge_strip.add(voxel_z_pair(v3, z3), delta_depth); // second edge
			} // for i
			if (num_ends > 0) {
				if (last_edge+2 < edge_strip.size()) edge_strip.resize(last_edge+2); // clip off extra points
				add_strip(edge_strip, 1, last_x, edge_y_start); // add if some edge elements
				edge_strip_len += edge_strip.size();
				++n_edge_strips;
			}
			add_strip(strip, 0, last_x, y_start);
			strip_len += strip.size();
			++n_strips;
		} // for n
	}
	cout << "strips: " << n_strips << " base, " << n_edge_strips << " edge" << endl;
	cout << "strip lengths: " << strip_len << " base, " << edge_strip_len << " edge" << endl;
	cout << "x_strip_map: " << x_strip_map.size() << endl;
	snow_draw.add_all_strips(snow_strips);
	snow_draw.finalize();
	snow_draw.show_stats();
}


bool snow_enabled() {

	return (temperature < W_FREEZE_POINT && !snow_draw.empty() && (display_mode & 0x02));
}


void gen_snow_coverage() {

	if (snow_depth <= 0.0 || num_snowflakes == 0) return; // disabled
	bool const vbo_supported(setup_gen_buffers());
	
	if (!vbo_supported) {
		cout << "Warning: VBOs not supported, so snow cannot be enabled." << endl;
		num_snowflakes = 0;
		return;
	}
	cout << "Determining Snow Coverage" << endl;
	voxel_map vmap;

	// setup voxel scales
	vox_delta.assign(VOXELS_PER_DIV/DX_VAL, VOXELS_PER_DIV/DY_VAL, 1.0/(max(DZ_VAL/VOXELS_PER_DIV, snow_depth)));
	
	if (read_snow_file) {
		RESET_TIME;
		vmap.read(snow_file);
		PRINT_TIME("Read Snow Voxel Map");
	}
	else {
		RESET_TIME;
		create_snow_map(vmap);
		PRINT_TIME("Build Snow Voxel Map");
	}
	if (write_snow_file) {
		RESET_TIME;
		vmap.write(snow_file);
		PRINT_TIME("Write Snow Voxel Map");
	}
	cout << "voxels: " << vmap.size() << endl;
	RESET_TIME;
	create_snow_strips(vmap);
	has_snow = snow_enabled();
	PRINT_TIME("Snow Strip Creation");
}


void reset_snow_vbos() {

	snow_draw.free_vbos();
}


void draw_snow() {

	has_snow = snow_enabled();
	if (!has_snow) return;
	static point last_lpos(all_zeros);
	point const lpos(get_light_pos());

	if (!shadow_map_enabled() && lpos != last_lpos) {
		RESET_TIME;
		update_cobj_tree();
		snow_draw.update_shadows();
		last_lpos = lpos;
		PRINT_TIME("Snow Shadow Calculation");
	}
	//RESET_TIME;
	shader_t s;

	if (!disable_shaders) {
		bool const use_smap(shadow_map_enabled());
		s.setup_enabled_lights();
		for (unsigned d = 0; d < 2; ++d) s.set_bool_prefix("no_normalize", !use_smap, d); // VS/FS
		s.set_prefix("#define USE_GOOD_SPECULAR", 1); // FS
		s.set_bool_prefix("use_shadow_map", use_smap, 1); // FS
		s.set_bool_prefix("use_texgen", 1, 0); // VS
		s.set_vert_shader("fog.part+texture_gen.part+per_pixel_lighting_textured");
		s.set_frag_shader("linear_fog.part+ads_lighting.part*+shadow_map.part*+per_pixel_lighting_textured");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		if (use_smap) set_smap_shader_for_all_lights(s);
	}
	set_specular(0.5, 50.0);
	set_color(SNOW_COLOR);
	plus_z.do_glNormal();
	point const camera(get_camera_pos());
	enable_blend();
	glDisable(GL_NORMALIZE);
	set_array_client_state(1, 0, 1, 0);
	select_texture(NOISE_TEX);
	setup_texgen(50.0, 50.0, 0.0, 0.0);
	snow_draw.draw();
	disable_textures_texgen();
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_NORMALIZE);
	disable_blend();
	set_specular(0.0, 1.0);
	s.end_shader();
	//PRINT_TIME("Snow Draw");
}


bool get_snow_height(point const &p, float radius, float &zval, vector3d &norm, bool crush_snow) {

	if (!has_snow || snow_strips.empty()) return 0;
	voxel_t const v(p);
	int const xval(v.p[0]), yval(v.p[1]);
	map<int, unsigned>::const_iterator it(x_strip_map.find(xval));
	if (it == x_strip_map.end()) return 0; // x-value not found
	assert(it->second < snow_strips.size());
	assert(snow_strips[it->second].xval == xval);
	
	for (unsigned i = it->second; i < snow_strips.size(); ++i) {
		strip_t const &s(snow_strips[i]);
		if (s.xval > xval)    break; // went too far in x
		assert(s.xval == xval);
		if (s.is_edge)        continue; // ignore edge strips for now
		if (s.y_start > yval) break; // went too far in y
		unsigned const len(s.get_size() >> 1); // divide by 2 because we want length in y
		assert(len >= 2);
		int const y_end(s.y_start + len - 1);
		if (yval > y_end)     continue; // not at y value yet
		unsigned const pos((yval - s.y_start) << 1);
		float const z(s.get_pos(pos).z);
		
		// Note: can get here for more than one strip if there is strip overlap in {y,z}
		if ((p.z - radius) < z && (p.z + radius) > z) {
			zval = z;
			norm = s.get_norm(pos);
			if (crush_snow) snow_draw.update_region(i, pos, s.get_size(), min(z, max((z - 0.25f*radius), (p.z - radius))));
			return 1;
		}
	}
	return 0;
}



