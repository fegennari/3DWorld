// 3D World - WaveFront Object File Reader
// by Frank Gennari
// 8/15/11
// Reference: http://en.wikipedia.org/wiki/Wavefront_.obj_file

#include "3DWorld.h"
#include <fstream>

using namespace std;


class object_file_reader {

protected:
	string filename;
	ifstream in;

	bool open_file() { // FIXME: what about multiple opens/reads?
		assert(!filename.empty());
		in.open(filename);
		if (in.good()) return 1;
		cerr << "Error: Could not open object file " << filename << endl;
		return 0;
	}

	void read_to_newline(ifstream &in_) const {
		// FIXME: what about '\' line wraps?
		in_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	void normalize_index(int &ix, unsigned vect_sz) const {
		assert(ix != 0);

		if (ix < 0) { // negative (relative) index
			assert(vect_sz >= (unsigned)(-ix));
			ix = vect_sz - ix;
		}
		else { // positive (absolute) index
			--ix; // specified starting from 1, but we want starting from 0
			assert((unsigned)ix < vect_sz);
		}
	}

public:
	object_file_reader(string const &fn) : filename(fn) {assert(!fn.empty());}

	bool read(vector<vector<point> > *ppts, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		vector<point> v; // vertices
		string s;

		while (in.good() && (in >> s)) {
			assert(!s.empty());

			if (s[0] == '#') { // comment
				read_to_newline(in); // ignore
			}
			else if (s == "v") { // vertex
				v.push_back(point());
			
				if (!(in >> v.back().x >> v.back().y >> v.back().z)) {
					cerr << "Error reading vertex from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "f") { // face
				if (ppts) ppts->push_back(vector<point>());
				int ix(0);

				while (in >> ix) { // read vertex index
					normalize_index(ix, v.size());
					if (ppts) ppts->back().push_back(v[ix]);

					if (in.get() == '/') {
						if (in >> ix) {} else in.clear(); // text coord index

						if (in.get() == '/') {
							if (in >> ix) {} else in.clear(); // normal index
						}
						else in.unget();
					}
					else in.unget();
				}
				in.clear();
			}
			else {
				read_to_newline(in); // ignore everything else
			}
		}
		PRINT_TIME("Polygons Load");
		if (verbose) cout << "v: " << v.size() << ", f: " << (ppts ? ppts->size() : 0) << endl;
		return 1;
	}
};


// ************************************************


colorRGB const def_color(0.0, 0.0, 0.0);

typedef map<string, unsigned> string_map_t;


struct polygon_t : public vector<vert_norm_tc> {

	void render(bool textured) const {
		assert(size() >= 3);

		for (const_iterator v = begin(); v != end(); ++v) {
			if (textured) glTexCoord2fv(v->t);
			v->n.do_glNormal();
			v->v.do_glVertex();
		}
	}
};


struct poly_group {

	vector<polygon_t> polygons;

	void add_poly(polygon_t const &poly) {polygons.push_back(poly);}
	void clear() {polygons.clear();}
	bool empty() const {return polygons.empty();}

	void render(bool textured) const {
		// FIXME: very inefficient - use glDrawArrays() or vbo with polygons split into triangles
		for (vector<polygon_t>::const_iterator p = polygons.begin(); p != polygons.end(); ++p) {
			glBegin(GL_POLYGON);
			p->render(textured);
			glEnd();
		}
	}
};


struct material_t {

	colorRGB ka, kd, ks, ke, tf;
	float ns, ni, alpha, tr;
	unsigned illum;
	int a_tid, d_tid, s_tid, alpha_tid, bump_tid;

	// geometry - does this go here or somewhere else?
	poly_group geom;

	material_t() : ka(def_color), kd(def_color), ks(def_color), ke(def_color), tf(def_color), ns(1.0), ni(1.0), alpha(1.0), tr(0.0),
		illum(2), a_tid(-1), d_tid(-1), s_tid(-1), alpha_tid(-1), bump_tid(-1) {}
	int get_render_texture() const {return d_tid;}

	void render(vector<texture_t> const &textures) const {
		if (geom.empty()) return; // nothing to do
		int const tex_id(get_render_texture());
		bool const textured(tex_id >= 0);

		if (textured) {
			assert((unsigned)tex_id < textures.size());
			assert(textures[tex_id].tid > 0);
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, textures[tex_id].tid);
		}
		if (alpha < 1.0 && ni != 1.0) {
			// set index of refraction (and reset it at the end)
		}
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		set_specular(spec_val, ns);
		set_color_a(colorRGBA(ka, alpha));
		set_color_d(colorRGBA(kd, alpha));
		set_color_e(colorRGBA(kd, alpha));
		geom.render(textured);
		set_color_e(BLACK);
		set_specular(0.0, 1.0);
		if (textured) glDisable(GL_TEXTURE_2D);
	}
};


class model3d {

	// geometry
	poly_group unbound_geom;

	// materials
	vector<material_t> materials;
	string_map_t mat_map; // maps material names to materials indexes
	set<string> undef_materials; // to reduce warning messages

	// textures
	vector<texture_t> textures;
	string_map_t tex_map; // maps texture filenames to texture indexes

public:
	unsigned num_materials(void) const {return materials.size();}

	material_t &get_material(int mat_id) {
		assert(mat_id >= 0 && (unsigned)mat_id < materials.size());
		return materials[mat_id];
	}

	// creation and query
	void add_polygon(polygon_t const &poly, int mat_id) {
		//assert(mat_id >= 0); // must be set/valid - FIXME: too strict?

		if (mat_id < 0) {
			unbound_geom.add_poly(poly);
		}
		else {
			assert((unsigned)mat_id < materials.size());
			materials[mat_id].geom.add_poly(poly);
		}
	}

	int get_material_ix(string const &material_name, string const &fn) {
		unsigned mat_id(0);
		string_map_t::const_iterator it(mat_map.find(material_name));

		if (it == mat_map.end()) {
			mat_id = materials.size();
			mat_map[material_name] = mat_id;
			materials.push_back(material_t());
		}
		else {
			cerr << "Warning: Redefinition of material " << material_name << " in file " << fn << endl;
			mat_id = it->second;
		}
		assert(mat_id < materials.size());
		return mat_id;
	}

	int find_material(string const &material_name) {
		string_map_t::const_iterator it(mat_map.find(material_name));

		if (it == mat_map.end()) {
			if (undef_materials.find(material_name) == undef_materials.end()) {
				cerr << "Error: Material " << material_name << " not found in any included material libraries" << endl;
				undef_materials.insert(material_name);
			}
			return -1; // return -1 on failure
		}
		assert(it->second < materials.size());
		return it->second;
	}

	unsigned create_texture(string const &fn, bool verbose) {
		string_map_t::const_iterator it(tex_map.find(fn));

		if (it != tex_map.end()) { // found (already loaded)
			assert(it->second < textures.size());
			return it->second;
		}
		unsigned const tid(textures.size());
		tex_map[fn] = tid;
		if (verbose) cout << "loading texture " << fn << endl;
		// type format width height wrap ncolors use_mipmaps name [bump_name [id [color]]]
		textures.push_back(texture_t(0, 4, 0, 0, 1, 3, 1, fn)); // always RGB targa wrapped+mipmap
		return tid; // can't fail
	}

	// clear/free
	void clear() {
		unbound_geom.clear();
		materials.clear();
		undef_materials.clear();
		mat_map.clear();
		free_textures();
		textures.clear();
		tex_map.clear();
	}

	void free_tids() {
		for (vector<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
			t->gl_delete();
		}
	}

	void free_textures() {
		for (vector<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
			t->free();
		}
	}

	// texture loading
	void ensure_texture_loaded(texture_t &t) const {
		if (!t.data) t.load(-1);
		assert(t.data);
	}

	void ensure_tid_loaded(int tid) {
		if (tid < 0) return; // not allocated
		assert((unsigned)tid < textures.size());
		ensure_texture_loaded(textures[tid]);
	}

	void ensure_tid_bound(int tid) {
		if (tid < 0) return; // not allocated
		assert((unsigned)tid < textures.size());
		textures[tid].check_init();
	}

	void load_all_used_tids() {
		for (vector<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			if (m->geom.empty()) continue;
			ensure_tid_loaded(m->get_render_texture()); // only one tid for now
			ensure_tid_loaded(m->alpha_tid);
			// FIXME: use alpha_tid
		}
	}

	void bind_all_used_tids() {
		load_all_used_tids();
		
		for (vector<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			if (m->geom.empty()) continue;
			ensure_tid_bound(m->get_render_texture()); // only one tid for now
		}
	}

	// rendering
	void render() { // const?
		bind_all_used_tids();
		unbound_geom.render(0);
		
		for (vector<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			m->render(textures);
		}
	}
};


// ************************************************


class object_file_reader_model : public object_file_reader {

	string rel_path;
	model3d &model;

	string open_include_file(string const &fn, string const &type, ifstream &in_inc) const {
		assert(!fn.empty());
		in_inc.open(fn); // try absolute path
		if (in_inc.good()) return fn;
		in_inc.clear();
		string const rel_fn(rel_path + fn);
		in_inc.open(rel_fn); // try relative path
		if (in_inc.good()) return rel_fn;
		cerr << "Error: Could not open " << type << " file " << fn << " or " << rel_fn << endl;
		return string();
	}

	string get_path(string const &fn) const {
		for (unsigned pos = fn.size(); pos > 0; --pos) {
			if (fn[pos-1] == '\\' || fn[pos-1] == '/') {
				return string(fn.begin(), fn.begin()+pos);
			}
		}
		return string();
	}

	int get_texture(string const &fn) {
		ifstream tex_in; // used only for determining file location
		string const fn_used(open_include_file(fn, "texture", tex_in));
		if (fn_used.empty()) return -1;
		tex_in.close();
		return model.create_texture(fn_used, 1);
	}

	void check_and_bind(int &tid, string const &tfn) {
		assert(tid < 0);
		tid = get_texture(tfn);
	}

public:
	object_file_reader_model(string const &fn, model3d &model_) : object_file_reader(fn), model(model_) {
		rel_path = get_path(fn);
	}

	bool load_mat_lib(string const &fn) { // FIXME: cache these files or reload every time?
		ifstream mat_in;
		if (open_include_file(fn, "material library", mat_in).empty()) return 0;
		cout << "loading material library " << fn << endl;
		int cur_mat_id(-1); // not set
		material_t *cur_mat(0);
		string s, tfn;

		while (mat_in.good() && (mat_in >> s)) {
			assert(!s.empty());

			if (s[0] == '#') { // comment - can we have comments here?
				read_to_newline(mat_in); // ignore
			}
			else if (s == "newmtl") { // new material
				string material_name;

				if (!(mat_in >> material_name)) {
					cerr << "Error reading material name" << endl;
					return 0;
				}
				cout << "Material " << material_name << endl; // FIXME: too verbose?
				cur_mat_id =  model.get_material_ix(material_name, fn);
				cur_mat    = &model.get_material(cur_mat_id);
			}
			else if (s == "Ka") {
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ka.R >> cur_mat->ka.G >> cur_mat->ka.B)) {cerr << "Error reading material Ka" << endl; return 0;}
			}
			else if (s == "Kd") {
				assert(cur_mat);
				if (!(mat_in >> cur_mat->kd.R >> cur_mat->kd.G >> cur_mat->kd.B)) {cerr << "Error reading material Kd" << endl; return 0;}
			}
			else if (s == "Ks") {
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ks.R >> cur_mat->ks.G >> cur_mat->ks.B)) {cerr << "Error reading material Ks" << endl; return 0;}
			}
			else if (s == "Ke") {
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ke.R >> cur_mat->ke.G >> cur_mat->ke.B)) {cerr << "Error reading material Ke" << endl; return 0;}
			}
			else if (s == "Ns") { // specular exponent
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ns)) {cerr << "Error reading material Ns" << endl; return 0;}
			}
			else if (s == "Ni") { // index of refraction
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ni)) {cerr << "Error reading material Ni" << endl; return 0;}
			}
			else if (s == "d") { // alpha
				assert(cur_mat);
				if (!(mat_in >> cur_mat->alpha)) {cerr << "Error reading material d" << endl; return 0;}
			}
			else if (s == "Tr") { // transmittance
				assert(cur_mat);
				if (!(mat_in >> cur_mat->tr)) {cerr << "Error reading material Tr" << endl; return 0;}
			}
			else if (s == "Tf") { // transmittion filter
				assert(cur_mat);
				if (!(mat_in >> cur_mat->tf.R >> cur_mat->tf.G >> cur_mat->tf.B)) {cerr << "Error reading material Tf" << endl; return 0;}
			}
			else if (s == "illum") { // 0 - 10
				assert(cur_mat);
				if (!(mat_in >> cur_mat->illum)) {cerr << "Error reading material Tr" << endl; return 0;}
			}
			else if (s == "map_Ka") {
				assert(cur_mat);
				if (!(mat_in >> tfn)) {cerr << "Error reading material map_Ka" << endl; return 0;}
				check_and_bind(cur_mat->a_tid, tfn);
			}
			else if (s == "map_Kd") {
				assert(cur_mat);
				if (!(mat_in >> tfn)) {cerr << "Error reading material map_Kd" << endl; return 0;}
				check_and_bind(cur_mat->d_tid, tfn);
			}
			else if (s == "map_Ks") {
				assert(cur_mat);
				if (!(mat_in >> tfn)) {cerr << "Error reading material map_Ks" << endl; return 0;}
				check_and_bind(cur_mat->s_tid, tfn);
			}
			else if (s == "map_d") {
				assert(cur_mat);
				if (!(mat_in >> tfn)) {cerr << "Error reading material map_d" << endl; return 0;}
				check_and_bind(cur_mat->alpha_tid, tfn);
			}
			else if (s == "map_bump" || s == "bump") { // should be ok if both are set
				assert(cur_mat);
				if (!(mat_in >> tfn)) {cerr << "Error reading material " << s << endl; return 0;}
				cur_mat->bump_tid = get_texture(tfn);
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in material library" << endl;
				return 0;
				//read_to_newline(mat_in); // ignore
			}
		}
		return 1;
	}

	bool read(vector<vector<point> > *ppts, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		int cur_mat_id(-1);
		vector<point> v; // vertices
		vector<vector3d> n; // normals
		vector<vector3d> tc; // texture coords
		string s, group_name, object_name, material_name, mat_lib, last_mat_lib;

		while (in.good() && (in >> s)) {
			assert(!s.empty());

			if (s[0] == '#') { // comment
				read_to_newline(in); // ignore
			}
			else if (s == "f") { // face
				if (ppts) ppts->push_back(vector<point>());
				polygon_t poly;
				int vix(0), tix(0), nix(0);
				bool has_norm(0);

				while (in >> vix) { // read vertex index
					normalize_index(vix, v.size());
					if (ppts) ppts->back().push_back(v[vix]);
					vector3d normal;
					vector3d tex_coord;

					if (in.get() == '/') {
						if (in >> tix) { // read text coord index
							normalize_index(tix, tc.size());
							tex_coord = tc[tix];
						}
						else {
							in.clear();
							tex_coord.x = tex_coord.y = 0.0; // better not be used
						}
						if (in.get() == '/') {
							if (in >> nix) { // read normal index
								normalize_index(nix, n.size());
								normal = n[nix];
								has_norm = 1;
							}
							else {
								in.clear();
								normal = zero_vector; // will be recalculated later
							}
						}
						else {
							in.unget();
						}
					}
					else {
						in.unget();
					}
					poly.push_back(vert_norm_tc(v[vix], normal, tex_coord.x, tex_coord.y));
				} // end while vertex
				if (!has_norm) { // calculate and set normal
					assert(poly.size() >= 3);
					vector3d const normal(cross_product((poly[1].v - poly[0].v), (poly[2].v - poly[0].v)).get_norm()); // backwards?
					for (unsigned i = 0; i < poly.size(); ++i) poly[i].n = normal;
				}
				in.clear();
				model.add_polygon(poly, cur_mat_id);
			}
			else if (s == "v") { // vertex
				v.push_back(point());
			
				if (!(in >> v.back().x >> v.back().y >> v.back().z)) {
					cerr << "Error reading vertex from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "vt") { // tex coord
				tc.push_back(vector3d());
			
				if (!(in >> tc.back().x >> tc.back().y)) {
					cerr << "Error reading texture coord from object file " << filename << endl;
					return 0;
				}
				if (!(in >> tc.back().z)) in.clear(); // optionally read a third tex coord
			}
			else if (s == "vn") { // normal
				n.push_back(vector3d());
			
				if (!(in >> n.back().x >> n.back().y >> n.back().z)) {
					cerr << "Error reading normal from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "l") { // line
				read_to_newline(in); // ignore
			}
			else if (s == "o") { // object definition
				if (!(in >> object_name)) {
					cerr << "Error reading object name from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "g") { // group
				if (!(in >> group_name)) {
					cerr << "Error reading group name from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "s") { // smoothing/shading (off/on or 0/1)
				read_to_newline(in); // ignore
			}
			else if (s == "usemtl") { // use material
				if (!(in >> material_name)) {
					cerr << "Error reading material from object file " << filename << endl;
					return 0;
				}
				cur_mat_id = model.find_material(material_name);
			}
			else if (s == "mtllib") { // material library
				if (!(in >> mat_lib)) {
					cerr << "Error reading material library from object file " << filename << endl;
					return 0;
				}
				if (mat_lib != last_mat_lib && !load_mat_lib(mat_lib)) {
					cerr << "Error reading material library file " << mat_lib << endl;
					return 0;
				}
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in object file " << filename << endl;
				return 0;
			}
		}
		PRINT_TIME("CXX Object Load");
		
		if (verbose) {
			cout << "v: " << v.size() << ", n: " << n.size() << ", tc: " << tc.size()
				 << ", f: " << (ppts ? ppts->size() : 0) << ", mat: " << model.num_materials() << endl;
		}
		return 1;
	}
};


bool read_object_file(char *filename, vector<vector<point> > &ppts, bool verbose) {

	if (1) {
		model3d model;
		object_file_reader_model reader(filename, model);
		return reader.read(&ppts, verbose);
	}
	else {
		object_file_reader reader(filename);
		return reader.read(&ppts, verbose);
	}
}


