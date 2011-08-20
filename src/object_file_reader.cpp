// 3D World - WaveFront Object File Reader
// by Frank Gennari
// 8/15/11
// Reference: http://en.wikipedia.org/wiki/Wavefront_.obj_file

#include "3DWorld.h"
#include "model3d.h"
#include <fstream>


extern bool recalc_model3d_normals;
extern model3ds all_models;


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

	bool read(vector<vector<point> > *ppts, geom_xform_t const &xf, bool verbose) {
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
				xf.xform_pos(v.back());
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
		return model.tm.create_texture(fn_used, 1);
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
			else if (s == "skip") { // skip this material
				assert(cur_mat);
				if (!(mat_in >> cur_mat->skip)) {cerr << "Error reading material skip" << endl; return 0;}
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in material library" << endl;
				return 0;
				//read_to_newline(mat_in); // ignore
			}
		}
		return 1;
	}


	struct vert_norm_tc_ix : public vert_norm_tc {
		unsigned ix;
		vert_norm_tc_ix(point const &v_, vector3d const &n_, float ts, float tt, unsigned ix_) : vert_norm_tc(v_, n_, ts, tt), ix(ix_) {}
	};

	struct poly_vix : public vector<vert_norm_tc_ix> {
		int mat_id;
		vector3d n;
		poly_vix(int mat_id_) : mat_id(mat_id_), n(zero_vector) {}
	};


	bool read(vector<vector<point> > *ppts, geom_xform_t const &xf, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		int cur_mat_id(-1);
		vector<point> v; // vertices
		vector<vector3d> n; // normals
		vector<vector3d> vn; // vertex normals
		vector<vector3d> tc; // texture coords
		deque<poly_vix> polys;
		string s, group_name, object_name, material_name, mat_lib, last_mat_lib;

		while (in.good() && (in >> s)) {
			assert(!s.empty());

			if (s[0] == '#') { // comment
				read_to_newline(in); // ignore
			}
			else if (s == "f") { // face
				polys.push_back(poly_vix(cur_mat_id));
				poly_vix &pv(polys.back());
				int vix(0), tix(0), nix(0);

				while (in >> vix) { // read vertex index
					normalize_index(vix, v.size());
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
					pv.push_back(vert_norm_tc_ix(v[vix], normal, tex_coord.x, tex_coord.y, vix));
				} // end while vertex
				in.clear();
				assert(pv.size() >= 3);
				vector3d const normal(cross_product((pv[1].v - pv[0].v), (pv[2].v - pv[0].v)).get_norm()); // backwards?
				pv.n = normal;

				for (unsigned i = 0; i < pv.size(); ++i) {
					assert((unsigned)pv[i].ix < vn.size());
					vn[pv[i].ix] += normal;
				}
			}
			else if (s == "v") { // vertex
				v.push_back(point());
				vn.push_back(zero_vector); // vertex normal
			
				if (!(in >> v.back().x >> v.back().y >> v.back().z)) {
					cerr << "Error reading vertex from object file " << filename << endl;
					return 0;
				}
				xf.xform_pos(v.back());
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
		PRINT_TIME("Object File Load");
		vntc_vect_t poly;

		for (vector<vector3d>::iterator i = vn.begin(); i != vn.end(); ++i) {
			i->normalize(); // average the normals
		}
		for (deque<poly_vix>::const_iterator i = polys.begin(); i != polys.end(); ++i) {
			poly.resize(i->size());

			for (unsigned j = 0; j < i->size(); ++j) {
				poly[j] = (*i)[j];
				if (!recalc_model3d_normals && poly[j].n != zero_vector) continue;
				assert((*i)[j].ix < vn.size());
				vector3d const &vert_norm(vn[(*i)[j].ix]);
				poly[j].n = ((fabs(dot_product(vert_norm, i->n)) < 0.75) ? i->n : vert_norm);
			}
			model.add_polygon(poly, i->mat_id, ppts);
		}
		PRINT_TIME("Model3d Build");

		if (verbose) {
			cout << "v: " << v.size() << ", n: " << n.size() << ", tc: " << tc.size()
				 << ", f: " << (ppts ? ppts->size() : 0) << ", mat: " << model.num_materials() << endl;
		}
		return 1;
	}
};


bool read_object_file(char *filename, vector<vector<point> > &ppts, geom_xform_t const &xf,
	int def_tid, colorRGBA const &def_c, bool load_models, bool verbose)
{
	if (load_models) {
		all_models.push_back(model3d(all_models.tm, def_tid, def_c));
		object_file_reader_model reader(filename, all_models.back());
		return reader.read(&ppts, xf, verbose);
	}
	else {
		object_file_reader reader(filename);
		return reader.read(&ppts, xf, verbose);
	}
}


