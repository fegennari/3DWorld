// 3D World - WaveFront Object File Reader
// by Frank Gennari
// 8/15/11
// Reference: http://en.wikipedia.org/wiki/Wavefront_.obj_file

#include "3DWorld.h"
#include "model3d.h"
#include "file_reader.h"
#include <stdint.h>
#include "fast_atof.h"


extern float model_auto_tc_scale;
extern model3ds all_models;

// hack to avoid slow multithreaded locking in getc()/ungetc() in MSVC++
#ifndef _getc_nolock
#define _getc_nolock   getc
#define _ungetc_nolock ungetc
#endif


bool base_file_reader::open_file(bool binary) {
	assert(!filename.empty());
	assert(!fp); // must call close_file() before reusing
	fp = fopen(filename.c_str(), (binary ? "rb" : "r"));
	if (!fp) {cerr << "Error: Could not open object file " << filename << endl;}
	return (fp != 0);
}

void base_file_reader::close_file() {
	if (fp) {fclose(fp);}
	fp = NULL;
}

void base_file_reader::unget_last_char(int c) {assert(fp); _ungetc_nolock(c, fp);}
int base_file_reader::get_char(FILE *fp_) const {return _getc_nolock(fp_);}

int base_file_reader::fast_atoi(char *str) const {
	//return atoi(str);
	assert(str && str[0] != 0);
	int v(0);
	bool is_neg(0);
	if (*str == '-') {is_neg = 1; ++str;} // negative

	for (; *str; ++str) {
		if (!fast_isdigit(*str)) return 0; // error
		v = 10*v + unsigned(*str - '0');
	}
	return (is_neg ? -v : v);
}

bool base_file_reader::read_int(int &v) {
	//return (fscanf(fp, "%i", &v) == 1);
	unsigned ix(0);
			
	while (1) {
		if (ix+1 >= MAX_CHARS) return 0; // buffer overrun
		char const c(get_next_char());
		if (ix == 0 && fast_isspace(c)) continue; // skip leading whitespace
		if (!fast_isdigit(c) && !(ix == 0 && c == '-')) {unget_last_char(c); break;} // non-integer character, unget it and finish
		buffer[ix++] = c;
	}
	if (ix == 0) return 0; // no integer characters were read
	buffer[ix] = 0; // add null terminator
	v = fast_atoi(buffer);
	return 1;
}

bool base_file_reader::read_string(char *s, unsigned max_len) {
	//return (fscanf(fp, "%s", s) == 1);
	unsigned ix(0);
			
	while (1) {
		if (ix+1 >= max_len) return 0; // buffer overrun
		char const c(get_next_char());
		if (fast_isspace(c)) {if (ix == 0) continue; else break;} // leading/trailing whitespace
		s[ix++] = c;
	}
	s[ix] = 0; // add null terminator
	return 1;
}


class object_file_reader : public base_file_reader {

protected:
	void normalize_index(int &ix, unsigned vect_sz) const {
		assert(ix != 0);

		if (ix < 0) { // negative (relative) index
			ix += vect_sz;
		}
		else { // positive (absolute) index
			--ix; // specified starting from 1, but we want starting from 0
		}
		assert((unsigned)ix < vect_sz);
	}

	bool read_point(point &p, unsigned req_num=3) {
		//return (fscanf(fp, "%f%f%f", &p.x, &p.y, &p.z) >= (int)req_num);
		for (unsigned i = 0; i < 3; ++i) {
			unsigned ix(0);
			
			while (1) {
				if (ix+1 >= MAX_CHARS) return 0; // buffer overrun
				char const c(get_next_char());
				if (fast_isspace(c)) {if (ix == 0) continue; else break;} // leading/trailing whitespace
				
				if (ix == 0 && !fast_isdigit(c) && c != '.' && c != '-') { // not a fp number
					unget_last_char(c);
					return ((i >= req_num) ? 1 : 0); // success if we read enough values
				}
				buffer[ix++] = c;
			}
			buffer[ix] = 0; // add null terminator
			p[i] = Assimp::fast_atof(buffer);
		}
		return 1;
	}

public:
	object_file_reader(string const &fn) : base_file_reader(fn) {}

	bool read(vector<coll_tquad> *ppts, geom_xform_t const &xf, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		cout << "Reading object file " << filename << endl;
		vector<point> v; // vertices
		char s[MAX_CHARS];
		polygon_t poly;

		while (read_string(s, MAX_CHARS)) {
			if (s[0] == '#') { // comment
				read_to_newline(fp); // ignore
			}
			else if (strcmp(s, "v") == 0) { // vertex
				v.push_back(point());
			
				if (!read_point(v.back())) {
					cerr << "Error reading vertex from object file " << filename << endl;
					return 0;
				}
				xf.xform_pos(v.back());
			}
			else if (strcmp(s, "f") == 0) { // face
				poly.resize(0);
				int ix(0);

				while (read_int(ix)) { // read vertex index
					normalize_index(ix, (unsigned)v.size());
					// only fill in the vertex (norm and tc will be unused)
					if (ppts) {poly.push_back(vert_norm_tc(v[ix], zero_vector, 0.0, 0.0));}
					int const c(get_next_char());

					if (c == '/') {
						read_int(ix); // text coord index, ok to fail
						int const c2(get_next_char());
						if (c2 == '/') {read_int(ix);} // normal index, ok to fail
						else {unget_last_char(c2);}
					}
					else {unget_last_char(c);}
				}
				if (ppts) {split_polygon(poly, *ppts, POLY_COPLANAR_THRESH);}
			}
			else {
				read_to_newline(fp); // ignore everything else
			}
		} // while
		PRINT_TIME("Polygons Load");
		if (verbose) cout << "v: " << v.size() << ", f: " << (ppts ? ppts->size() : 0) << endl;
		return 1;
	}
};


// ************************************************


string model_from_file_t::open_include_file(string const &fn, string const &type, ifstream &in_inc) const {
	assert(!fn.empty());
	// try absolute path
	in_inc.open(fn);
	if (in_inc.good()) return fn;
	// try relative path
	in_inc.clear();
	string const rel_fn(rel_path + fn);
	in_inc.open(rel_fn);
	if (in_inc.good()) return rel_fn;
	// try texture directory
	in_inc.clear();
	string const tex_fn("textures/" + fn);
	in_inc.open(tex_fn);
	if (in_inc.good()) return tex_fn;
	cerr << "Error: Could not open " << type << " file " << fn << ", " << rel_fn << ", or " << tex_fn << endl;
	return string();
}

string model_from_file_t::get_path(string const &fn) const {
	for (unsigned pos = (unsigned)fn.size(); pos > 0; --pos) {
		if (fn[pos-1] == '\\' || fn[pos-1] == '/') {return string(fn.begin(), fn.begin()+pos);}
	}
	return string();
}

int model_from_file_t::get_texture(string const &fn, bool is_alpha_mask, bool verbose, bool invert_alpha, bool wrap, bool mirror) {
	int const existing_tid(texture_lookup(fn));
	if (existing_tid >= 0) {return (BUILTIN_TID_START + existing_tid);}
	ifstream tex_in; // used only for determining file location
	string const fn_used(open_include_file(fn, "texture", tex_in));
	if (fn_used.empty()) return -1;
	tex_in.close();
	return model.tmgr.create_texture(fn_used, is_alpha_mask, verbose, invert_alpha, wrap, mirror);
}



class object_file_reader_model : public object_file_reader, public model_from_file_t {

	bool read_map_name(ifstream &in, string &name) const {
		if (!(in >> name))  {return 0;}
		assert(!name.empty());
		if (name[0] == '-') {return 1;} // option, return and let the parser deal with it
		read_to_newline(in, &name);
		return 1;
	}

	bool read_color_rgb(ifstream &in, string const &name, colorRGB &color) const {
		// Note: does not support "xyz" CIE-XYZ color space specification
		if (!(in >> color.R)) {cerr << "Error reading material " << name << endl; return 0;}
		if (!(in >> color.G >> color.B)) { // G and B are optional
			in.clear(); // clear error bits
			color.G = color.B = color.R; // grayscale
		}
		return 1;
	}

public:
	object_file_reader_model(string const &fn, model3d &model_) : object_file_reader(fn), model_from_file_t(fn, model_) {}

	bool load_mat_lib(string const &fn) { // Note: could cache filename, but seems to never be included more than once
		ifstream mat_in;
		if (open_include_file(fn, "material library", mat_in).empty()) return 0;
		cout << "loading material library " << fn << endl;
		int cur_mat_id(-1); // not set
		material_t *cur_mat(0);
		string s, tfn;

		while (mat_in.good() && (mat_in >> s)) {
			assert(!s.empty());
			transform(s.begin(), s.end(), s.begin(), tolower); // convert all characters to lowercase

			if (s[0] == '#') { // comment - can we have comments here?
				read_to_newline(mat_in); // ignore
			}
			else if (s == "newmtl") { // new material
				string material_name;

				if (!(mat_in >> material_name)) {
					cerr << "Error reading material name" << endl;
					return 0;
				}
				if (verbose) cout << "Material " << material_name << endl; // maybe too verbose?
				cur_mat_id =  model.get_material_ix(material_name, fn);
				cur_mat    = &model.get_material(cur_mat_id);
			}
			else if (s == "ka") {
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ka", cur_mat->ka)) return 0;
			}
			else if (s == "kd") {
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Kd", cur_mat->kd)) return 0;
			}
			else if (s == "ks") {
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ks", cur_mat->ks)) return 0;
			}
			else if (s == "ke") {
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ke", cur_mat->ke)) return 0;
			}
			else if (s == "ns") { // specular exponent
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ns)) {cerr << "Error reading material Ns" << endl; return 0;}
			}
			else if (s == "ni") { // index of refraction
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ni)) {cerr << "Error reading material Ni" << endl; return 0;}
			}
			else if (s == "d") { // dissolve, treated as alpha
				assert(cur_mat);
				if (!(mat_in >> cur_mat->alpha)) {cerr << "Error reading material d" << endl; return 0;}
			}
			else if (s == "tr") { // transmittance
				assert(cur_mat);
				if (!(mat_in >> cur_mat->tr)) {cerr << "Error reading material Tr" << endl; return 0;}
			}
			else if (s == "tf") { // transmittion filter
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Tf", cur_mat->tf)) return 0;
			}
			else if (s == "illum") { // 0 - 10
				assert(cur_mat);
				if (!(mat_in >> cur_mat->illum)) {cerr << "Error reading material Tr" << endl; return 0;}
			}
			else if (s == "sharpness") { // Note: unused
				assert(cur_mat);
				float sharpness; // Note: unused
				if (!(mat_in >> sharpness)) {cerr << "Error reading material sharpness" << endl; return 0;}
			}
			// Note: may want to support the -clamp option to clamp textures
			else if (s == "map_ka") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ka" << endl; return 0;}
				check_and_bind(cur_mat->a_tid, tfn, 0, verbose);
			}
			else if (s == "map_kd") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Kd" << endl; return 0;}
				check_and_bind(cur_mat->d_tid, tfn, 0, verbose); // invert=0, wrap=1, mirror=0
			}
			else if (s == "map_ks") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ks" << endl; return 0;}
				check_and_bind(cur_mat->s_tid, tfn, 0, verbose);
			}
			else if (s == "map_d") { // dissolve
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_d" << endl; return 0;}
				check_and_bind(cur_mat->alpha_tid, tfn, 1, verbose);
			}
			else if (s == "map_bump" || s == "bump") { // should be ok if both are set
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material " << s << endl; return 0;}
				
				if (tfn == "-bm") { // this is the only material sub-option we handle, and the scale is ignored
					float scale(1.0);
					if (!(mat_in >> scale) || !read_map_name(mat_in, tfn)) {cerr << "Error reading material " << s << " with -bm" << endl; return 0;}
				}
				cur_mat->bump_tid = get_texture(tfn, 0, verbose); // can be set from both map_bump and bump
			}
			else if (s == "map_refl") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_refl" << endl; return 0;}
				check_and_bind(cur_mat->refl_tid, tfn, 0, verbose);
			}
			else if (s == "map_Ns") { // Note: unused
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ns" << endl; return 0;}
			}
			//else if (s == "map_aat") {} // toggle antialiasing
			//else if (s == "decal") {} // modifies color
			//else if (s == "dist") {} // displacement
			else if (s == "skip") { // skip this material
				assert(cur_mat);
				if (!(mat_in >> cur_mat->skip)) {cerr << "Error reading material skip" << endl; return 0;}
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in material library. Skipping line." << endl;
				read_to_newline(mat_in); // ignore
				//return 0;
			}
		} // while
		return 1;
	}


	bool load_from_model3d_file(bool verbose) {
		RESET_TIME;

		if (!model.read_from_disk(filename)) {
			cerr << "Error reading model3d file " << filename << endl;
			return 0;
		}
		PRINT_TIME("Model3d File Load");
		set<string> mat_lib_fns;
		model.get_all_mat_lib_fns(mat_lib_fns);
		
		for (set<string>::const_iterator i = mat_lib_fns.begin(); i != mat_lib_fns.end(); ++i) {
			if (!load_mat_lib(*i)) {
				cerr << "Error reading material library file " << *i << endl;
				//return 0;
			}
		}
		model.load_all_used_tids();

		if (verbose) {
			cout << "bcube: " << model.get_bcube().str() << endl;
			cout << "model stats: "; model.show_stats();
		}
		PRINT_TIME("Model3d Load");
		return 1;
	}


	bool read(geom_xform_t const &xf, bool recalc_normals, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		cout << "Reading object file " << filename << endl;
		unsigned const block_size = (1 << 18); // 256K
		int cur_mat_id(-1);
		unsigned smoothing_group(0), prev_smoothing_group(0), num_faces(0), num_objects(0), num_groups(0), obj_group_id(0);
		vector<point> v; // vertices
		vector<vector3d> n; // normals
		vector<counted_normal> vn; // vertex normals
		vector<point2d<float> > tc; // texture coords
		deque<poly_data_block> pblocks;
		set<string> loaded_mat_libs;
		char s[MAX_CHARS];
		string material_name, mat_lib, group_name, object_name;
		tc.push_back(point2d<float>(0.0, 0.0)); // default tex coords
		n.push_back(zero_vector); // default normal
		unsigned approx_line(0);

		while (read_string(s, MAX_CHARS)) {
			++approx_line;

			if (s[0] == 0) {
				cout << "empty/unparseable line?" << endl;
				continue;
			}
			else if (s[0] == '#') { // comment
				read_to_newline(fp); // ignore
			}
			else if (strcmp(s, "f") == 0) { // face
				model.mark_mat_as_used(cur_mat_id);

				if (pblocks.empty() || pblocks.back().pts.size() >= block_size || smoothing_group != prev_smoothing_group) { // create a new block
					if (!pblocks.empty()) {
						remove_excess_cap(pblocks.back().polys);
						remove_excess_cap(pblocks.back().pts);
					}
					pblocks.push_back(poly_data_block());
					prev_smoothing_group = smoothing_group;
				}
				poly_data_block &pb(pblocks.back());
				pb.polys.push_back(poly_header_t(cur_mat_id, obj_group_id));
				unsigned &npts(pb.polys.back().npts);
				unsigned const pix((unsigned)pb.pts.size());
				int vix(0), tix(0), nix(0);

				while (read_int(vix)) { // read vertex index
					normalize_index(vix, (unsigned)v.size());
					vntc_ix_t vntc_ix(vix, 0, 0);
					int const c(get_next_char());

					if (c == '/') {
						if (read_int(tix)) { // read text coord index
							normalize_index(tix, (unsigned)tc.size()-1); // account for tc[0]
							vntc_ix.tix = tix+1; // account for tc[0]
						}
						int const c2(get_next_char());

						if (c2 == '/') {
							if (read_int(nix) && !recalc_normals) { // read normal index
								normalize_index(nix, (unsigned)n.size()-1); // account for n[0]
								vntc_ix.nix = nix+1; // account for n[0]
							} // else the normal will be recalculated later
						}
						else {unget_last_char(c2);}
					}
					else {unget_last_char(c);}
					pb.pts.push_back(vntc_ix);
					++npts;
				} // end while vertex
				assert(npts >= 3);
				vector3d &normal(pb.polys.back().n);
				
				for (unsigned i = pix; i < pix+npts-2; ++i) { // find a nonzero normal
					normal = cross_product((v[pb.pts[i+1].vix] - v[pb.pts[i].vix]), (v[pb.pts[i+2].vix] - v[pb.pts[i].vix])); // backwards?
					// if we disable this normalize() we will weight normal contributions by polygon area,
					// but we have to change the code below and it causes problems with vertex uniquing
					normal.normalize();
					if (normal != zero_vector) break; // got a good normal
				}
				if (recalc_normals) {
					for (unsigned i = pix; i < pix+npts; ++i) {
						unsigned const vix(pb.pts[i].vix);
						assert((unsigned)vix < vn.size());

						if (vn[vix].is_valid() && dot_product(normal, vn[vix].get_norm()) < 0.25) { // normals in disagreement
							vn[vix] = zero_vector; // zero it out so that it becomes invalid later
						}
						else {
							vn[vix].add_normal(normal);
						}
					}
				}
			}
			else if (strcmp(s, "v") == 0) { // vertex
				v.push_back(point());
				if (recalc_normals) vn.push_back(counted_normal()); // vertex normal
			
				if (!read_point(v.back())) {
					cerr << "Error reading vertex from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				xf.xform_pos(v.back());
			}
			else if (strcmp(s, "vt") == 0) { // tex coord
				point tc3d;
			
				if (!read_point(tc3d, 2)) {
					cerr << "Error reading texture coord from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				tc.push_back(point2d<float>(tc3d.x, tc3d.y)); // discard tc3d.z
			}
			else if (strcmp(s, "vn") == 0) { // normal
				vector3d normal;
			
				if (!read_point(normal)) {
					cerr << "Error reading normal from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				if (!recalc_normals) {
					xf.xform_pos_rm(normal);
					n.push_back(normal);
				}
			}
			else if (strcmp(s, "l") == 0) { // line
				read_to_newline(fp); // ignore
			}
			else if (strcmp(s, "o") == 0) { // object definition
				read_str_to_newline(fp, object_name); // can be empty?
				++num_objects;
				++obj_group_id;
			}
			else if (strcmp(s, "g") == 0) { // group
				read_str_to_newline(fp, group_name); // can be empty
				++num_groups;
				++obj_group_id;
			}
			else if (strcmp(s, "s") == 0) { // smoothing/shading (off/on or 0/1)
				if (fscanf(fp, "%u", &smoothing_group) != 1) {
					if (!read_string(s, MAX_CHARS) || strcmp(s, "off") != 0) {
						cerr << "Error reading smoothing group from object file " << filename << " near line " << approx_line << endl;
						return 0;
					}
					smoothing_group = 0;
				}
			}
			else if (strcmp(s, "usemtl") == 0) { // use material
				read_str_to_newline(fp, material_name);

				if (material_name.empty()) {
					cerr << "Error reading material from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				cur_mat_id = model.find_material(material_name);
			}
			else if (strcmp(s, "mtllib") == 0) { // material library
				read_str_to_newline(fp, mat_lib);

				if (mat_lib.empty()) {
					cerr << "Error reading material library from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				if (loaded_mat_libs.find(mat_lib) == loaded_mat_libs.end()) { // mtllib not yet loaded
					if (!load_mat_lib(mat_lib)) { // can materials be redefined with different mat_libs?
						cerr << "Error reading material library file " << mat_lib << " near line " << approx_line << endl;
						//return 0;
					}
					loaded_mat_libs.insert(mat_lib);
				}
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in object file " << filename << " near line " << approx_line << endl;
				read_to_newline(fp); // ignore this line
				//return 0;
			}
		} // while
		remove_excess_cap(v);
		remove_excess_cap(n);
		remove_excess_cap(tc);
		remove_excess_cap(vn);
		PRINT_TIME("Object File Load");
		model.load_all_used_tids(); // need to load the textures here to get the colors
		PRINT_TIME("Model Texture Load");
		size_t const num_blocks(pblocks.size());
		model3d::proc_counted_normals(vn); // if recalc_normals

		while (!pblocks.empty()) {
			poly_data_block const &pd(pblocks.back());
			unsigned pix(0);
			polygon_t poly;
			vntc_map_t vmap[2]; // {triangles, quads}
			vntct_map_t vmap_tan[2]; // {triangles, quads}

			for (vector<poly_header_t>::const_iterator j = pd.polys.begin(); j != pd.polys.end(); ++j) {
				poly.resize(j->npts);
				
				for (unsigned p = 0; p < j->npts; ++p) {
					vntc_ix_t const &V(pd.pts[pix+p]);
					vector3d normal;

					if (recalc_normals) {
						assert(V.vix < vn.size());
						normal = ((j->n != zero_vector && !vn[V.vix].is_valid()) ? j->n : vn[V.vix]);
					}
					else {
						assert(V.nix < n.size());
						normal = n[V.nix];
						if (normal == zero_vector) normal = j->n;
					}
					assert(V.vix < v.size() && V.tix < tc.size());
					point2d<float> tcoord;

					if (V.tix == 0 && model_auto_tc_scale > 0.0) { // generate tc since it wasn't read from the file
						unsigned const dim(get_max_dim(normal)), dimx((dim == 0) ? 1 : 0), dimy((dim == 2) ? 1 : 2); // looks better for brick textures on walls
						tcoord.x = model_auto_tc_scale*v[V.vix][dimx];
						tcoord.y = model_auto_tc_scale*v[V.vix][dimy];
					}
					else {tcoord = tc[V.tix];}
					poly[p] = vert_norm_tc(v[V.vix], normal, tcoord.x, tcoord.y);
				}
				num_faces += model.add_polygon(poly, vmap, vmap_tan, j->mat_id, j->obj_id);
				pix += j->npts;
			}
			pblocks.pop_back();
		}
		model.optimize(); // optimize vertices and remove excess capacity
		PRINT_TIME("Model3d Build");
		
		if (verbose) {
			size_t const nn(recalc_normals ? vn.size() : n.size());
			cout << "verts: " << v.size() << ", normals: " << nn << ", tcs: " << tc.size() << ", faces: " << num_faces << ", objects: " << num_objects
				 << ", groups: " << num_groups << ", blocks: " << num_blocks << endl;
			cout << "bcube: " << model.get_bcube().str() << endl << "model stats: "; model.show_stats();
		}
		return 1;
	}
};


void check_obj_file_ext(string const &filename, string const &ext) {
	if (ext != "obj") {cout << "Warning: Attempting to read file '" << filename << "' with extension '" << ext << "' as an object file." << endl;}
}


bool write_model3d_file(string const &base_fn, model3d &cur_model) {

	RESET_TIME;
	assert(base_fn.size() > 4);
	string out_fn(base_fn.begin(), base_fn.end()-4); // strip off the '.obj'
	out_fn += ".model3d";
	cur_model.bind_all_used_tids(); // need to force tangent vector calculation
				
	if (!cur_model.write_to_disk(out_fn)) {
		cerr << "Error writing model3d file " << out_fn << endl;
		return 0;
	}
	PRINT_TIME("Model3d Write");
	return 1;
}


bool read_3ds_file_model(string const &filename, model3d &model, geom_xform_t const &xf, bool use_vertex_normals, bool verbose);
bool read_3ds_file_pts(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, colorRGBA const &def_c, bool verbose);


bool read_model_file(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, int def_tid, colorRGBA const &def_c,
	int reflective, float metalness, bool load_models, bool recalc_normals, bool write_file, bool verbose)
{
	string const ext(get_file_extension(filename, 0, 1));
	std::locale::global(std::locale("C"));
	setlocale(LC_ALL, "C");

	if (load_models) {
		all_models.push_back(model3d(all_models.tmgr, def_tid, def_c, reflective, metalness));
		model3d &cur_model(all_models.back());

		if (ext == "3ds") {
			if (!read_3ds_file_model(filename, cur_model, xf, recalc_normals, verbose)) return 0; // recalc_normals is always true
		}
		else { // object/model3d file
			object_file_reader_model reader(filename, cur_model);

			if (ext == "model3d") {
				//assert(xf == geom_xform_t()); // xf is ignored, assumed to be already applied; use transforms with loaded model3d files
				if (!reader.load_from_model3d_file(verbose)) return 0;
			}
			else {
				check_obj_file_ext(filename, ext);
				if (!reader.read(xf, recalc_normals, verbose)) return 0;
				if (write_file && !write_model3d_file(filename, cur_model)) return 0;
			}
		}
		if (ppts) {get_cur_model_polygons(*ppts);}
		return 1;
	}
	else {
		if (ext == "3ds") {
			return read_3ds_file_pts(filename, ppts, xf, def_c, verbose);
		}
		else {
			check_obj_file_ext(filename, ext);
			object_file_reader reader(filename);
			return reader.read(ppts, xf, verbose);
		}
	}
}


