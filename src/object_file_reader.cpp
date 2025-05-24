// 3D World - WaveFront Object File Reader
// by Frank Gennari
// 8/15/11
// Reference: http://en.wikipedia.org/wiki/Wavefront_.obj_file

#include "3DWorld.h"
#include "model3d.h"
#include "file_reader.h"
#include <stdint.h>
#include <algorithm> // for transform()
#include <cctype> // for tolower()
#include "fast_atof.h"
#include "format_text.h"


extern bool use_obj_file_bump_grayscale, model_calc_tan_vect, enable_model_animations, enable_spec_map, enable_shine_map;
extern float model_auto_tc_scale, model_mat_lod_thresh;
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
	if (fp) {checked_fclose(fp);}
	fp = NULL;
}

bool base_file_reader::read_int(int &v) {
	bool first_char(1), is_neg(0);
	v = 0;

	while (1) {
		char const c(get_char(fp));
		if (first_char && fast_isspace(c)) continue; // skip leading whitespace
		if (first_char && c == '-') {is_neg = 1; first_char = 0; continue;} // negative
		if (!fast_isdigit(c)) {unget_last_char(c); break;} // non-integer character, unget it and finish
		v = 10*v + int(c - '0');
		first_char = 0;
	}
	if (first_char) return 0; // no integer characters were read
	if (is_neg) {v = -v;}
	return 1;
}
bool base_file_reader::read_uint(unsigned &v) {
	int temp(-1);
	if (!read_int(temp) || temp < 0) return 0;
	v = temp; // cast to unsigned
	return 1;
}

bool base_file_reader::read_string(char *s, unsigned max_len) {
	unsigned ix(0);
			
	while (1) {
		if (ix+1 >= max_len) return 0; // buffer overrun
		char const c(get_char(fp));
		if (c == EOF) break;
		if (fast_isspace(c)) {
			if (ix == 0) continue; // leading whitespace
			if (c == '\n') {unget_last_char(c);} // preserve the newline
			break; // trailing whitespace
		}
		s[ix++] = c;
	}
	if (ix == 0) return 0; // nothing was read
	s[ix] = 0; // add null terminator
	return 1;
}


class object_file_reader : public base_file_reader {

	bool invalid_index_warned;

protected:
	void handle_invalid_zero_ref_index(int &ix) {
		if (ix == -1) {
			if (!invalid_index_warned) {
				cerr << "Error: Invalid zero index in object file" << endl;
				//assert(0);
				invalid_index_warned = 1;
			}
			++ix;
		}
	}
	void normalize_index(int &ix, unsigned vect_sz) {
		int const input_ix(ix);
		if (ix < 0) { // negative (relative) index
			ix += vect_sz;
		}
		else { // positive (absolute) index
			--ix; // specified starting from 1, but we want starting from 0
		}
		handle_invalid_zero_ref_index(ix);
		if ((unsigned)ix >= vect_sz) {cout << TXT(input_ix) << TXT(ix) << TXT(vect_sz) << endl;}
		assert((unsigned)ix < vect_sz);
	}
	bool read_float(float &val) {
		unsigned ix(0);

		while (1) {
			if (ix+1 >= MAX_CHARS) return 0; // buffer overrun
			char const c(get_char(fp));
			if (c == EOF) break;
			if (fast_isspace(c)) {if (ix == 0) continue; else break;} // leading/trailing whitespace

			if (ix == 0 && !fast_isdigit(c) && c != '.' && c != '-') { // not a fp number
				unget_last_char(c);
				return 0;
			}
			buffer[ix++] = c;
		}
		buffer[ix] = 0; // add null terminator
		val = Assimp::fast_atof(buffer);
		return 1;
	}
	bool read_point(point &p, unsigned req_num=3) {
		//return (fscanf(fp, "%f%f%f", &p.x, &p.y, &p.z) >= (int)req_num);
		for (unsigned i = 0; i < 3; ++i) {
			if (!read_float(p[i])) {return ((i >= req_num) ? 1 : 0);} // success if we read enough values
		}
		return 1;
	}
	int read_optional_color_RGB(colorRGB &c) { // return value: 0=no color read, 1=color read, 2=error
		float val(0.0);
		if (!read_float(val)) return 0; // no more numbers to read
		c.R = val;
		return ((read_float(c.G) && read_float(c.B)) ? 1 : 2); // success or error
	}

public:
	object_file_reader(string const &fn) : base_file_reader(fn), invalid_index_warned(0) {}

	bool read(vector<coll_tquad> *ppts, geom_xform_t const &xf, bool verbose) {
		RESET_TIME;
		if (!open_file(1)) return 0; // binary mode is faster
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
					if (ppts) {poly.emplace_back(v[ix], zero_vector, 0.0, 0.0);}
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

string model_from_file_t::get_path(string const &fn) {
	for (unsigned pos = (unsigned)fn.size(); pos > 0; --pos) {
		if (fn[pos-1] == '\\' || fn[pos-1] == '/') {return string(fn.begin(), fn.begin()+pos);}
	}
	return string();
}

int model_from_file_t::get_texture(string const &fn, bool is_alpha_mask, bool verbose, bool invert_alpha, bool wrap, bool mirror, bool force_grayscale) {
	int const existing_tid(texture_lookup(fn));
	if (existing_tid >= 0) {return (BUILTIN_TID_START + existing_tid);}
	ifstream tex_in; // used only for determining file location
	string const fn_used(open_include_file(fn, "texture", tex_in));
	if (fn_used.empty()) return -1;
	tex_in.close();
	return model.tmgr.create_texture(fn_used, is_alpha_mask, verbose, invert_alpha, wrap, mirror, force_grayscale);
}

void model_from_file_t::check_and_bind(int &tid, string const &tfn, bool is_alpha_mask, bool verbose, bool invert_alpha, bool wrap, bool mirror) {
	if (tid >= 0) {
		cerr << "Warning: Duplicate specification of object file material texture " << tfn << " ignored" << endl;
		return;
	}
	tid = get_texture(tfn, is_alpha_mask, verbose, invert_alpha, wrap, mirror, 0);
}


bool startswith(string const &str, string const &prefix) {
	if (prefix.size() > str.size()) return 0;
	return std::equal(prefix.begin(), prefix.end(), str.begin());
}
bool endswith(string const &str, string const &suffix) {
	if (suffix.size() > str.size()) return 0;
	return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

class object_file_reader_model : public object_file_reader, public model_from_file_t {

	bool had_empty_mat_error;

	bool read_map_name(ifstream &in, string &name, float *scale=nullptr) {
		if (!(in >> name)) {return 0;} // no name read (EOF?)
		assert(!name.empty());

		if (name[0] == '-') {
			if (name == "-bm") { // this is the only material sub-option we handle, and the scale is ignored
				float scale_(1.0);
				if (!(in >> scale_)) return 0; // read the scale
				if (scale != nullptr) {*scale = scale_;}
				if (!(in >> name)) {return 0;} // now read the actual name
			}
			else {
				cerr << "Map option " << name << " is not supported." << endl;
				return 0;
			}
		}
		read_to_newline(in, &name);
		// Note: sometimes the -bm option can be after the filename instead
		// however, this is difficult to differentiate from a filename with whitespace and hyphen, so we don't support it and require the user to update their file
		// unless it's "-bm 1", which seems to be the common case
		if (endswith(name, " -bm 1")) {name = name.substr(0, name.size()-6);}
		return 1;
	}
	bool read_color_rgb(ifstream &in, string const &name, colorRGB &color, string const &material_name) const {
		// Note: does not support "xyz" CIE-XYZ color space specification
		if (!(in >> color.R)) {cerr << "Error reading material " << name << " for " << material_name << endl; return 0;}
		if (!(in >> color.G >> color.B)) { // G and B are optional
			in.clear(); // clear error bits
			color.G = color.B = color.R; // grayscale
		}
		//color.set_valid_color(); // colors should be in the [0.0, 1.0] range?
		return 1;
	}

public:
	object_file_reader_model(string const &fn, model3d &model_) : object_file_reader(fn), model_from_file_t(fn, model_), had_empty_mat_error(0) {}

	bool load_mat_lib(string const &fn) { // Note: could cache filename, but seems to never be included more than once
		ifstream mat_in;
		if (open_include_file(fn, "material library", mat_in).empty()) return 0;
		cout << "Loading material library " << fn << endl;
		int cur_mat_id(-1); // not set
		material_t *cur_mat(0);
		string s, tfn, material_name;

		while (mat_in.good() && (mat_in >> s)) {
			assert(!s.empty());
			transform(s.begin(), s.end(), s.begin(), ::tolower); // convert all characters to lowercase

			if (s[0] == '#') { // comment - can we have comments here?
				read_to_newline(mat_in); // ignore
			}
			else if (s == "newmtl") { // new material
				material_name.clear();
				read_to_newline(mat_in, &material_name);

				if (material_name.empty()) {
					cerr << "Error reading material name" << endl;
					return 0;
				}
				if (verbose) {cout << "Material " << material_name << endl;} // maybe too verbose?
				cur_mat_id =  model.get_material_ix(material_name, fn);
				assert(cur_mat_id >= 0); // must be valid
				cur_mat    = &model.get_material(cur_mat_id);
			}
			else if (s == "ka") { // ambient color
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ka", cur_mat->ka, material_name)) return 0;
			}
			else if (s == "kd") { // diffuse color
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Kd", cur_mat->kd, material_name)) return 0;
			}
			else if (s == "ks") { // specular color
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ks", cur_mat->ks, material_name)) return 0;
			}
			else if (s == "ke") { // emissive color
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Ke", cur_mat->ke, material_name)) return 0;
			}
			else if (s == "ns") { // specular exponent
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ns)) {cerr << "Error reading material Ns for " << material_name << endl; return 0;}
			}
			else if (s == "ni") { // index of refraction
				assert(cur_mat);
				if (!(mat_in >> cur_mat->ni)) {cerr << "Error reading material Ni for " << material_name << endl; return 0;}
			}
			else if (s == "d") { // dissolve, treated as alpha
				assert(cur_mat);
				if (!(mat_in >> cur_mat->alpha)) {cerr << "Error reading material d for " << material_name << endl; return 0;}
			}
			else if (s == "tr") { // transmittance
				assert(cur_mat);
				if (!(mat_in >> cur_mat->tr)) {cerr << "Error reading material Tr for " << material_name << endl; return 0;}
			}
			else if (s == "tf") { // transmittion filter - stored but unused
				assert(cur_mat);
				if (!read_color_rgb(mat_in, "Tf", cur_mat->tf, material_name)) return 0;
			}
			else if (s == "km") { // bump strengh - ignored
				assert(cur_mat);
				float km(0.0);
				if (!(mat_in >> km)) {cerr << "Error reading material km for " << material_name << endl; return 0;}
			}
			else if (s == "illum") { // 0 - 10
				assert(cur_mat);
				if (!(mat_in >> cur_mat->illum)) {cerr << "Error reading material Tr for " << material_name << endl; return 0;}
			}
			else if (s == "sharpness") { // Note: unused
				assert(cur_mat);
				float sharpness; // Note: unused
				if (!(mat_in >> sharpness)) {cerr << "Error reading material sharpness for " << material_name << endl; return 0;}
			}
			else if (s == "clamp") { // Note: unused, normally is On/Off
				assert(cur_mat);
				string val; // Note: unused
				if (!(mat_in >> val)) {cerr << "Error reading material clamp for " << material_name << endl; return 0;}
				// Note: may want to support the -clamp option to clamp textures
			}
			else if (s == "map_ka") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ka for " << material_name << endl; return 0;}
				check_and_bind(cur_mat->a_tid, tfn, 0, verbose);
			}
			else if (s == "map_kd") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Kd for " << material_name << endl; return 0;}
				check_and_bind(cur_mat->d_tid, tfn, 0, verbose); // invert=0, wrap=1, mirror=0
			}
			else if (s == "map_ks") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ks for " << material_name << endl; return 0;}
				if (enable_spec_map) {check_and_bind(cur_mat->s_tid, tfn, 0, verbose);}
			}
			else if (s == "map_ns") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_Ns for " << material_name << endl; return 0;}
				if (enable_shine_map) {check_and_bind(cur_mat->ns_tid, tfn, 0, verbose);}
			}
			else if (s == "map_d") { // dissolve
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_d for " << material_name << endl; return 0;}
				check_and_bind(cur_mat->alpha_tid, tfn, 1, verbose);
			}
			else if (s == "map_bump" || s == "bump" || s == "norm") { // should be ok if more than one are set; 3DWorld auto detects grayscale bump maps vs. RGB normal maps
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material " << s << " for " << material_name << endl; return 0;}
				cur_mat->bump_tid = get_texture(tfn, 0, verbose, 0, 1, 0, use_obj_file_bump_grayscale); // can be set from map_bump, bump, and norm
			}
			else if (s == "map_refl") {
				assert(cur_mat);
				if (!read_map_name(mat_in, tfn)) {cerr << "Error reading material map_refl for " << material_name << endl; return 0;}
				check_and_bind(cur_mat->refl_tid, tfn, 0, verbose);
			}
			// unsupported
			else if (s == "kt") {unhandled(s, mat_in);} // transmission color
			else if (s == "map_aat") {unhandled(s, mat_in);} // toggle antialiasing
			else if (s == "decal"  ) {unhandled(s, mat_in);} // modifies color
			else if (s == "disp"   ) {unhandled(s, mat_in);} // displacement
			// PBR parameters, mostly unsupported
			else if (s == "metalness" || s == "pm") {
				// Note: metalness has been added for 3DWorld and is not in the original obj file format;
				assert(cur_mat);
				if (!(mat_in >> cur_mat->metalness)) {cerr << "Error reading material metalness for " << material_name << endl; return 0;}
			}
			else if (s == "pr"    ) {unhandled(s, mat_in);} // roughness
			else if (s == "ps"    ) {unhandled(s, mat_in);} // sheen
			else if (s == "pc"    ) {unhandled(s, mat_in);} // clearcoat thickness
			else if (s == "pcr"   ) {unhandled(s, mat_in);} // clearcoat roughness
			else if (s == "aniso" ) {unhandled(s, mat_in);} // anisotropy
			else if (s == "anisor") {unhandled(s, mat_in);} // anisotropy rotation
			else if (s == "refl"  ) {unhandled(s, mat_in);} // reflection? seen in mailbox OBJ file
			// PBR maps
			else if (s == "map_pr") {unhandled(s, mat_in);} // metallic
			else if (s == "map_pm") {unhandled(s, mat_in);} // roughness
			else if (s == "map_ps") {unhandled(s, mat_in);} // sheen
			else if (s == "map_ke") {unhandled(s, mat_in);} // emissive
			// 3DWorld extensions
			else if (s == "skip") { // skip this material
				assert(cur_mat);
				if (!(mat_in >> cur_mat->skip)) {cerr << "Error reading material skip for " << material_name << endl; return 0;}
			}
			else {
				cerr << "Error: Undefined entry '" << s << "' in material library. Skipping line." << endl;
				read_to_newline(mat_in); // ignore
				//return 0;
			}
		} // while
		return 1;
	}
	void unhandled(string const &s, ifstream &mat_in) {
		cerr << "Warning: Unhandled entry '" << s << "' in material library. Skipping line." << endl;
		read_to_newline(mat_in); // ignore
	}

	bool load_from_model3d_file(bool verbose) {
		if (!model.read_from_disk(filename)) {
			cerr << format_red("Error reading model3d file " + filename) << endl;
			return 0;
		}
		{ // open a scope
			timer_t timer("Model3d File Load");
			set<string> mat_lib_fns;
			model.get_all_mat_lib_fns(mat_lib_fns);
		
			for (auto i = mat_lib_fns.begin(); i != mat_lib_fns.end(); ++i) {
				if (!load_mat_lib(*i)) {
					cerr << "Error reading material library file " << *i << endl;
					//return 0;
				}
			}
		}
		model.load_all_used_tids(); // optional
		if (verbose) {model.show_stats();}
		return 1;
	}


	bool try_load_mat_lib(string const &mat_lib, set<string> &loaded_mat_libs, unsigned approx_line) {
		if (loaded_mat_libs.find(mat_lib) == loaded_mat_libs.end()) { // mtllib not yet loaded
			if (!load_mat_lib(mat_lib)) { // can materials be redefined with different mat_libs?
				istringstream iss(mat_lib);
				bool ret(1);
				string str;
				
				while (iss >> str) { // try to split by whitespace, in case there were multiple mtllib filenames on the same line (instead of treating as whitespace in the filename)
					if (str == mat_lib) { // no whitespace
						cerr << "Error reading material library file " << str << " near line " << approx_line << endl;
						ret = 0;
					}
					else if (!try_load_mat_lib(str, loaded_mat_libs, approx_line)) {ret = 0;}
				}
				loaded_mat_libs.insert(mat_lib); // mark as loaded, even if load failed
				return ret;
			}
			loaded_mat_libs.insert(mat_lib);
		}
		return 1;
	}

	bool read(geom_xform_t const &xf, int recalc_normals, bool verbose) {
		if (!open_file(1)) return 0; // binary mode is faster
		cout << "Reading object file " << filename << endl;
		RESET_TIME;
		unsigned const block_size = (1 << 18); // 256K
		int cur_mat_id(-1);
		unsigned smoothing_group(0), prev_smoothing_group(0), num_faces(0), num_objects(0), num_groups(0), obj_group_id(0);
		vector<point> v; // vertices
		vector<vector3d> n; // normals
		// weighted_normal can also be used, but doesn't work well; see face_weight_avg mode selected by recalc_normals==2
		vector<counted_normal> vn; // vertex normals
		vector<point2d<float> > tc; // texture coords
		vector<colorRGB> colors; // vertex colors
		deque<poly_data_block> pblocks;
		set<string> loaded_mat_libs;
		char s[MAX_CHARS];
		string material_name, mat_lib, group_name, object_name;
		tc.push_back(point2d<float>(0.0, 0.0)); // default tex coords
		n.push_back(zero_vector); // default normal
		unsigned approx_line(0);
		bool is_textured(0), had_npts_error(0);

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
				unsigned const pix((unsigned)pb.pts.size()), pts_start(pb.pts.size());
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
				if (npts < 3) {
					if (!had_npts_error) {cerr << "Error near line " << approx_line << ": face has only " << npts << " vertices." << endl; had_npts_error = 1;}
					pb.pts.resize(pts_start);
					pb.polys.pop_back(); // remove pts and polygon
					continue; // skip it
				}
				vector3d &normal(pb.polys.back().n);
				
				for (unsigned i = pix; i < pix+npts-2; ++i) { // find a nonzero normal
					normal = cross_product((v[pb.pts[i+1].vix] - v[pb.pts[i].vix]), (v[pb.pts[i+2].vix] - v[pb.pts[i].vix])); // backwards?
					// if we disable this normalize() we will weight normal contributions by polygon area,
					// but we have to change the code below and it causes problems with vertex uniquing
					normal.normalize();
					if (normal != zero_vector) break; // got a good normal
				}
				if (recalc_normals) {
					bool const face_weight_avg(recalc_normals == 2 && (npts == 3 || npts == 4)); // only works for quads and triangles
					float face_area(0.0);

					if (face_weight_avg) {
						point face_pts[4];
						for (unsigned i = 0; i < npts; ++i) {face_pts[i] = v[pb.pts[i+pix].vix];}
						face_area = polygon_area(face_pts, npts);
					}
					for (unsigned i = pix; i < pix+npts; ++i) {
						unsigned const vix(pb.pts[i].vix);
						assert((unsigned)vix < vn.size());
						bool const using_texgen(is_textured && model_auto_tc_scale > 0.0 && pb.pts[i].tix == 0);

						if (vn[vix].is_valid() && (using_texgen || dot_product(normal, vn[vix].get_norm()) < 0.25)) { // normals in disagreement (or using texgen)
							vn[vix] = zero_vector; // zero it out so that it becomes invalid later
						}
						else if (face_weight_avg) {vn[vix].add_normal(face_area*normal);} // face weighted average
						else {vn[vix].add_normal(normal);} // unweighted average of normals
					}
				}
			}
			else if (strcmp(s, "v") == 0) { // vertex
				v.push_back(point());
				if (recalc_normals) {vn.push_back(counted_normal());} // vertex normal
			
				if (!read_point(v.back())) {
					cerr << "Error reading vertex from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				colorRGB color;
				int const color_ret(read_optional_color_RGB(color));
				if (color_ret == 2) {cerr << "Error reading vertex color from object file " << filename << " near line " << approx_line << endl; return 0;}
				else if (color_ret == 1) {
					if (colors.empty()) {colors.resize(v.size()-1, WHITE);} // pad colors up to this point with white
					colors.push_back(color);
				}
				else if (!colors.empty()) {colors.push_back(WHITE);} // color not specified, and in colors mode, pad with white
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
				if (!read_uint(smoothing_group)) {
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
					if (!had_empty_mat_error) {cerr << "Error reading material from object file " << filename << " near line " << approx_line << endl;}
					had_empty_mat_error = 1;
					return 0;
				}
				cur_mat_id = model.find_material(material_name);
				
				if (cur_mat_id >= 0) { // material was valid
					int const tid(model.get_material(cur_mat_id).d_tid);
					is_textured = (tid >= 0 && model.tmgr.get_tex_avg_color(tid) != WHITE); // no texture, or all white texture
				}
			}
			else if (strcmp(s, "mtllib") == 0) { // material library
				read_str_to_newline(fp, mat_lib);

				if (mat_lib.empty()) {
					cerr << "Error reading material library from object file " << filename << " near line " << approx_line << endl;
					return 0;
				}
				if (!try_load_mat_lib(mat_lib, loaded_mat_libs, approx_line)) {
					//return 0; // nonfatal
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
		remove_excess_cap(colors);
		PRINT_TIME("Object File Load");
		model.load_all_used_tids(); // need to load the textures here to get the colors
		size_t const num_blocks(pblocks.size());
		model3d::proc_model_normals(vn, recalc_normals); // if recalc_normals

		while (!pblocks.empty()) {
			poly_data_block const &pd(pblocks.back());
			unsigned pix(0);
			polygon_t poly;
			vntc_map_t vmap[2]; // {triangles, quads}
			vntct_map_t vmap_tan[2]; // {triangles, quads}
			colorRGBA color;

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
					if (!colors.empty()) {assert(V.vix < colors.size()); color += colors[V.vix];}
				} // for p
				if (!colors.empty()) {color = color/j->npts; color.A = 1.0;} // uses average vertex color for each face/polygon, with alpha=1.0
				// TODO: use color; model3d doesn't support per-vertex colors
				num_faces += model.add_polygon(poly, vmap, vmap_tan, j->mat_id, j->obj_id);
				pix += j->npts;
			} // for j
			pblocks.pop_back();
		}
		model.finalize(); // optimize vertices, remove excess capacity, compute bounding cube, subdivide, generate LOD blocks
		PRINT_TIME("Model3d Build");
		
		if (verbose) {
			size_t const nn(recalc_normals ? vn.size() : n.size());
			cout << "verts: " << v.size() << ", normals: " << nn << ", tcs: " << tc.size() << ", colors: " << colors.size() << ", faces: " << num_faces
				 << ", objects: " << num_objects << ", groups: " << num_groups << ", blocks: " << num_blocks << endl;
			model.show_stats();
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
	if (model_calc_tan_vect) {cur_model.calc_tangent_vectors();} // tangent vectors are needed for writing
				
	if (!cur_model.write_to_disk(out_fn)) {
		cerr << "Error writing model3d file " << out_fn << endl;
		return 0;
	}
	PRINT_TIME("Model3d Write");
	return 1;
}


bool read_3ds_file_model(string const &filename, model3d &model, geom_xform_t const &xf, int use_vertex_normals, bool verbose);
bool read_3ds_file_pts(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, colorRGBA const &def_c, bool verbose);
bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, string const &anim_name, int recalc_normals, bool verbose);

bool const ALWAYS_USE_ASSIMP = 0;

// recalc_normals: 0=no, 1=yes, 2=face_weight_avg
bool load_model_file(string const &filename, model3ds &models, geom_xform_t const &xf, string const &anim_name, int def_tid, colorRGBA const &def_c,
	int reflective, float metalness, float lod_scale, int recalc_normals, int group_cobjs_level, bool write_file, bool verbose, uint64_t rev_winding_mask)
{
	if (filename.empty()) return 0; // can't be loaded
	string const ext(get_file_extension(filename, 0, 1));
	models.push_back(model3d(filename, models.tmgr, def_tid, def_c, reflective, metalness, lod_scale, recalc_normals, group_cobjs_level));
	model3d &cur_model(models.back());

	if (!ALWAYS_USE_ASSIMP && ext == "3ds") {
		if (!read_3ds_file_model(filename, cur_model, xf, recalc_normals, verbose)) {models.pop_back(); return 0;} // recalc_normals is always true
		//if (write_file && !write_model3d_file(filename, cur_model)) return 0; // Note: doesn't work because there's no mtllib file
	}
	else if (ext == "model3d") {
		//assert(xf == geom_xform_t()); // xf is ignored, assumed to be already applied; use transforms with loaded model3d files
		if (!object_file_reader_model(filename, cur_model).load_from_model3d_file(verbose)) {models.pop_back(); return 0;}
	}
	else if (!ALWAYS_USE_ASSIMP && ext == "obj") {
		check_obj_file_ext(filename, ext);
		//test_other_obj_loader(filename); // placeholder for testing other object file loaders (tinyobjloader, assimp, etc.)
		if (!object_file_reader_model(filename, cur_model).read(xf, recalc_normals, verbose)) {models.pop_back(); return 0;}
		if (write_file && !write_model3d_file(filename, cur_model)) return 0; // don't need to pop the model
	}
	else { // not a built-in supported format, try using assimp if compiled in
		if (!read_assimp_model(filename, cur_model, xf, anim_name, recalc_normals, verbose)) {models.pop_back(); return 0;}
	}
	if (model_mat_lod_thresh > 0.0) {cur_model.compute_area_per_tri();} // used for TT LOD/distance culling
	cur_model.reverse_winding_order(rev_winding_mask);
	return 1;
}

// Note: assimp reader not supported in this flow
bool read_model_file(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, int def_tid, colorRGBA const &def_c,
	int reflective, float metalness, float lod_scale, bool load_models, int recalc_normals, int group_cobjs_level, bool write_file, bool verbose)
{
	setlocale(LC_ALL, "C"); // optimization for obj file reading?

	if (load_models) {
		// anim_name is used for debugging and general model loading; loaders that expect custom named animations don't call this function
		string const anim_name(enable_model_animations ? "default" : "");
		if (!load_model_file(filename, all_models, xf, anim_name, def_tid, def_c, reflective, metalness,
			lod_scale, recalc_normals, group_cobjs_level, write_file, verbose)) return 0;
		if (ppts) {get_cur_model_polygons(*ppts);}
		return 1;
	}
	else {
		string const ext(get_file_extension(filename, 0, 1));
		if (ext == "3ds") {return read_3ds_file_pts(filename, ppts, xf, def_c, verbose);}
		else {
			check_obj_file_ext(filename, ext);
			object_file_reader reader(filename);
			return reader.read(ppts, xf, verbose);
		}
	}
}


