// 3D World - WaveFront Object File Reader
// by Frank Gennari
// 8/15/11

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


typedef vector<vert_norm_tc> polygon_t;


struct poly_group {
	vector<polygon_t> polygons;
};


struct material_t {

	texture_t texture;
	// FIXME: write

	// geometry - does this go here or somewhere else?
	poly_group geom;
};


class model3d {

	poly_group unbound_geom;
	vector<material_t> materials;
	typedef map<string, unsigned> mat_map_t;
	mat_map_t mat_map; // maps material names to materials
	set<string> undef_materials; // to reduce warning messages

public:
	unsigned num_materials(void) const {return materials.size();}

	material_t &get_material(int mat_id) {
		assert(mat_id >= 0 && (unsigned)mat_id < materials.size());
		return materials[mat_id];
	}

	void add_polygon(polygon_t const &poly, int mat_id) {
		//assert(mat_id >= 0); // must be set/valid - FIXME: too strict?

		if (mat_id < 0) {
			unbound_geom.polygons.push_back(poly);
		}
		else {
			assert((unsigned)mat_id < materials.size());
			materials[mat_id].geom.polygons.push_back(poly);
		}
	}

	int get_material_ix(string const &material_name, string const &fn) {
		unsigned mat_id(0);
		mat_map_t::const_iterator it(mat_map.find(material_name));

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
		mat_map_t::const_iterator it(mat_map.find(material_name));

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

	bool load_texture(string const &fn, unsigned mat_id) {
		assert(mat_id < materials.size());
		material_t &mat(materials[mat_id]);
		assert(!mat.texture.data); // should not yet be loaded - FIXME: too strict?
		mat.texture.free();
		// type format width height wrap ncolors use_mipmaps name [bump_name [id [color]]]
		mat.texture = texture_t(0, 4, 0, 0, 1, 3, 1, fn); // always RGB targa wrapped+mipmap
		mat.texture.load(-1);
		return 1; // can't fail
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

	bool load_texture(string const &fn, unsigned mat_id) {
		ifstream tex_in; // used only for determining file location
		string const fn_used(open_include_file(fn, "texture", tex_in));
		if (fn_used.empty()) return 0;
		in.close();
		cout << "loading texture " << fn_used << endl; // FIXME: too verbose?
		return model.load_texture(fn_used, mat_id);
	}

	bool load_mat_lib(string const &fn) { // FIXME: cache these files or reload every time?
		ifstream mat_in;
		if (open_include_file(fn, "material library", mat_in).empty()) return 0;
		cout << "loading material library " << fn << endl;
		int cur_mat_id(-1); // not set
		string s;

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
				cur_mat_id = model.get_material_ix(material_name, fn);
			}
			/*else if (s == "Ka") {

			}
			else if (s == "Kd") {

			}
			else if (s == "Ks") {

			}
			else if (s == "Ke") {

			}
			else if (s == "Ns") {

			}
			else if (s == "Ni") {

			}
			else if (s == "d") {

			}
			else if (s == "Tr") {

			}
			else if (s == "Tf") {

			}
			else if (s == "illum") {

			}
			else if (s == "map_Ka") {

			}
			else if (s == "map_Kd") {

			}
			else if (s == "map_d") {

			}
			else if (s == "map_bump") {

			}
			else if (s == "bump") {

			}*/
			else {
				//cerr << "Error: Undefined entry '" << s << "' in material library" << endl;
				//return 0;
				read_to_newline(mat_in); // ignore
			}
		}
		return 1;
	}

public:
	object_file_reader_model(string const &fn, model3d &model_) : object_file_reader(fn), model(model_) {
		rel_path = get_path(fn);
	}

	bool read(vector<vector<point> > *ppts, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		int cur_mat_id(-1);
		vector<point> v; // vertices
		vector<vector3d> n; // normals
		vector<vector3d> tc; // texture coords
		string s, group_name, material_name, mat_lib, last_mat_lib;

		while (in.good() && (in >> s)) {
			assert(!s.empty());

			if (s[0] == '#') { // comment
				read_to_newline(in); // ignore
			}
			else if (s == "f") { // face
				if (ppts) ppts->push_back(vector<point>());
				polygon_t poly;
				int vix(0), tix(0), nix(0);

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
							}
							else {
								in.clear();
								// calc normal
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
				read_to_newline(in); // ignore
			}
			else if (s == "g") { // group
				if (!(in >> group_name)) {
					cerr << "Error reading group name from object file " << filename << endl;
					return 0;
				}
			}
			else if (s == "s") { // smoothing
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


