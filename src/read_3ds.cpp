// 3D World - 3D Studio Max File Reader
// by Frank Gennari
// 11/1/2014
// Reference: http://www.spacesimulator.net/wiki/index.php?title=Tutorials:3ds_Loader

#include "3DWorld.h"
#include "model3d.h"
#include "file_reader.h"


class file_reader_3ds : public base_file_reader {

protected:
	geom_xform_t cur_xf;
	float master_scale;
	string name; // unused?

	struct face_t {
		unsigned short ix[3], flags;
		int mat;
		face_t() : flags(0), mat(-1) {}
	};

	long get_end_pos(unsigned read_len) {return (ftell(fp) + read_len - 6);}

	bool read_data(void *dest, size_t sz, size_t count, char const *const str) {
		size_t const nread(fread(dest, sz, count, fp));
		if (nread == count) return 1;
		if (str != nullptr) {cerr << "Error reading 3DS file " << str << " data: expected " << count << " elements of size " << sz << " but got " << nread << endl;}
		return 0;
	}

	bool read_chunk_header(unsigned short &chunk_id, unsigned &chunk_len, bool top_level=0) {
		// Note: we pass null at the top level to suppress the error message during EOF checking
		if (!read_data(&chunk_id, 2, 1, (top_level ? nullptr : "chunk id"))) return 0;
		if (verbose) {printf("ChunkID: %x\n", chunk_id);}
		if (!read_data(&chunk_len, 4, 1, "chunk length")) return 0;
		assert(chunk_len < (1<<31)); // sanity check
		if (verbose) {printf("ChunkLength: %u\n", chunk_len);}
		return 1;
	}

	template<typename T> void ensure_size(vector<T> &v, unsigned num) {
		if (v.empty()) {v.resize(num);} else {assert(v.size() == num);}
	}

	bool read_vertex_block(vector<vert_tc_t> &verts) {
		unsigned short num;
		if (!read_data(&num, sizeof(unsigned short), 1, "vertex size")) return 0;
		ensure_size(verts, num);
		if (verbose) {printf("Number of vertices: %u\n", num);}

		for (int i = 0; i < num; i++) {
			if (!read_data(&verts[i].v.x, sizeof(float), 3, "vertex xyz")) return 0;
			verts[i].v *= master_scale;
			cur_xf.xform_pos(verts[i].v);
		}
		return 1;
	}

	bool read_mapping_block(vector<vert_tc_t> &verts) {
		unsigned short num;
		if (!read_data(&num, sizeof(unsigned short), 1, "mapping size")) return 0;
		ensure_size(verts, num);
		if (verbose) {printf("Number of tex coords: %u\n", num);}

		for (int i = 0; i < num; i++) {
			if (!read_data(verts[i].t, sizeof(float), 2, "mapping uv")) return 0;
		}
		return 1;
	}

	bool read_faces(vector<face_t> &faces) {
		unsigned short num;
		if (!read_data(&num, sizeof(unsigned short), 1, "number of faces")) return 0;
		ensure_size(faces, num);
		if (verbose) {printf("Number of polygons: %u\n", num);}

		for (unsigned i = 0; i < num; i++) {
			// flags: 3 LSB are edge orderings (AB BC CA): ignored
			if (!read_data(faces[i].ix, sizeof(unsigned short), 4, "face data")) return 0;
			std::swap(faces[i].ix[0], faces[i].ix[2]); // reverse triangle vertex ordering to agree with 3DWorld coordinate system
			//cout << "flags: " << faces[i].flags << endl;
		}
		return 1;
	}

	bool read_null_term_string(string &str) {
		while (1) {
			char l_char;
			if (!read_data(&l_char, 1, 1, "string char")) return 0;
			if (l_char == '\0') return 1;
			str.push_back(l_char); // add the null terminator or not?
		}
		return 0; // never gets here
	}

	bool read_color(colorRGB &color) {
		unsigned short chunk_id;
		unsigned chunk_len;
		if (!read_chunk_header(chunk_id, chunk_len)) return 0;

		if (chunk_id == 0x0010) {
			assert(chunk_len == 18);
			return read_data(&color.R, sizeof(float), 3, "RGB float color");
		}
		else if (chunk_id == 0x0011) {
			assert(chunk_len == 9);
			color_wrapper cw;
			if (!read_data(cw.c, sizeof(unsigned char), 3, "RGB 24-bit color")) return 0;
			color = cw.get_c3();
			return 1;
		}
		return 0; // invalid chunk_id
	}

	bool read_percentage(unsigned read_len, float &val) {
		unsigned short chunk_id, ival;
		unsigned chunk_len;
		long const end_pos(get_end_pos(read_len));
		if (!read_chunk_header(chunk_id, chunk_len)) return 0;
		
		switch (chunk_id) {
		case 0x0030: // short percentage
			if (!read_data(&ival, sizeof(unsigned short), 1, "short percentage")) return 0;
			val = ival/65535.0;
			break;
		case 0x0031: // float percentage
			if (!read_data(&val, sizeof(float), 1, "float percentage")) return 0;
			break;
		default:
			assert(0);
		} // end switch
		assert(ftell(fp) == end_pos);
		return 1;
	}

	void skip_chunk(unsigned chunk_len) {fseek(fp, chunk_len-6, SEEK_CUR);}

	// return value: 0 = error, 1 = processed, 2 = can't handle/skip
	virtual int proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {return 2;}

public:
	file_reader_3ds(string const &fn) : base_file_reader(fn), master_scale(1.0) {}
	virtual ~file_reader_3ds() {}

	bool read(geom_xform_t const &xf, bool verbose_) {
		RESET_TIME;
		verbose = verbose_;
		cur_xf = xf;
		if (!open_file(1)) return 0; // binary file
		cout << "Reading 3DS file " << filename << endl;
		unsigned short chunk_id;
		unsigned chunk_len;

		while (read_chunk_header(chunk_id, chunk_len, 1)) { // read each chunk from the file 
			switch (chunk_id) {
				// MAIN3DS: Main chunk, contains all the other chunks; length: 0 + sub chunks
			case 0x4d4d:
				break;
				// EDIT3DS: 3D Editor chunk, objects layout info; length: 0 + sub chunks
			case 0x3d3d:
				break;
				// EDIT_OBJECT: Object block, info for each object; length: len(object name) + sub chunks
			case 0x4000:
				if (!read_null_term_string(name)) return 0;
				break;
				// master scale factor
			case 0x0100:
				if (!read_data(&master_scale, sizeof(float), 1, "master scale")) return 0;
				break;
			default: // send to derived class reader, and skip chunk if it isn't handled
				{
					int const ret(proc_other_chunks(chunk_id, chunk_len));
					if (ret == 0) return 0; // error
					else if (ret == 2) {skip_chunk(chunk_len);} // skip
				}
			} // end switch
		}
		return 1;
	}
};


// ************************************************


class file_reader_3ds_triangles : public file_reader_3ds {

	colorRGBA def_color;
	vector<coll_tquad> *ppts;

	virtual int proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {
		assert(ppts != nullptr);

		switch (chunk_id) {
			// OBJ_TRIMESH: Triangular mesh, contains chunks for 3d mesh info; length: 0 + sub chunks
		case 0x4100:
			return read_mesh(chunk_len); // handled
		} // end switch
		return 2; // skip
	}

	bool read_mesh(unsigned read_len) {
		unsigned short chunk_id;
		unsigned chunk_len;
		long const end_pos(get_end_pos(read_len));
		vector<vert_tc_t> verts;
		vector<face_t> faces;

		while (ftell(fp) < end_pos) { // read each chunk from the file 
			if (!read_chunk_header(chunk_id, chunk_len)) return 0;

			switch (chunk_id) {
				// TRI_FACEL1: Polygons (faces) list
				// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
			case 0x4120:
				if (!read_faces(faces)) return 0;
				break;
				// TRI_VERTEXL: Vertices list
				// Chunk Length: 1 x unsigned short (# vertices) + 3 x float (vertex coordinates) x (# vertices) + sub chunks
			case 0x4110:
				if (!read_vertex_block(verts)) return 0;
				break;
				// TRI_MAPPINGCOORS: Vertices list
				// Chunk Length: 1 x unsigned short (# mapping points) + 2 x float (mapping coordinates) x (# mapping points) + sub chunks
			case 0x4140:
				if (!read_mapping_block(verts)) return 0;
				break;
			default:
				skip_chunk(chunk_len);
			} // end switch
		} // end while
		assert(ftell(fp) == end_pos);
		triangle tri;

		for (vector<face_t>::const_iterator i = faces.begin(); i != faces.end(); ++i) {
			for (unsigned n = 0; n < 3; ++n) {
				unsigned const ix(i->ix[n]);
				assert(ix < verts.size());
				tri.pts[n] = verts[ix].v;
			}
			ppts->push_back(coll_tquad(tri, def_color));
		}
		return 1;
	}

public:
	file_reader_3ds_triangles(string const &fn) : file_reader_3ds(fn), ppts(nullptr) {}

	bool read(vector<coll_tquad> *ppts_, geom_xform_t const &xf, colorRGBA const &def_c, bool verbose) {
		assert(ppts_ != nullptr);
		ppts = ppts_;
		def_color = def_c;
		return file_reader_3ds::read(xf, verbose);
	}
};


// ************************************************


class file_reader_3ds_model : public file_reader_3ds, public model_from_file_t {

	bool use_vertex_normals;
	unsigned obj_id;

	virtual int proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {

		// TODO: smoothing groups and shininess
		switch (chunk_id) {
			// OBJ_TRIMESH: Triangular mesh, contains chunks for 3d mesh info; length: 0 + sub chunks
		case 0x4100:
			return read_mesh(chunk_len); // handled

		case 0xAFFF: // material
			{
				// since the material properties may be defined before its name, we can't get the material by name and fill it in;
				// instead, we create a temporary material, fill it in, then look it up by name and overwrite the material in the model with cur_mat
				material_t cur_mat("", filename);
				if (!read_material(chunk_len, cur_mat)) {cerr << "Error reading material " << cur_mat.name << endl; return 0;}
				int const cur_mat_id(model.get_material_ix(cur_mat.name, filename, 0)); // FIXME: okay_if_exists???
				model.get_material(cur_mat_id) = cur_mat;
				return 1; // handled
			}
		} // end switch
		return 2; // skip
	}

	bool read_mesh(unsigned read_len) {
		unsigned short chunk_id;
		unsigned chunk_len;
		long const end_pos(get_end_pos(read_len));
		vector<vert_tc_t> verts;
		vector<face_t> faces;
		typedef map<int, vector<unsigned short> > face_mat_map_t;
		face_mat_map_t face_materials;

		while (ftell(fp) < end_pos) { // read each chunk from the file 
			if (!read_chunk_header(chunk_id, chunk_len)) return 0;

			switch (chunk_id) {
				// TRI_FACEL1: Polygons (faces) list
				// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
			case 0x4120:
				if (!read_faces(faces)) return 0;
				break;
			// faces data
			case 0x4130: // Faces Material: asciiz name, short nfaces, short face_ids
				{
					// read and process material name
					string mat_name;
					if (!read_null_term_string(mat_name)) return 0;
					int const mat_id(model.get_material_ix(mat_name, filename, 1));
					model.mark_mat_as_used(mat_id);
					vector<unsigned short> &faces_mat(face_materials[mat_id]);
					assert(faces_mat.empty());
					// read and process face materials
					unsigned short num;
					if (!read_data(&num, sizeof(unsigned short), 1, "number of faces for material")) return 0;
					assert(num > 0); // too strong?
					faces_mat.resize(num);
					if (!read_data(&faces_mat.front(), sizeof(unsigned short), num, "faces for material")) return 0;
					if (verbose) {cout << "Material " << mat_name << " is used for " << num << " faces" << endl;}
					break;
				}
			case 0x4150: // Smoothing Group List
				// nfaces*4bytes: Long int where the nth bit indicates if the face belongs to the nth smoothing group
				skip_chunk(chunk_len); // FIXME
				break;
				// TRI_VERTEXL: Vertices list
				// Chunk Length: 1 x unsigned short (# vertices) + 3 x float (vertex coordinates) x (# vertices) + sub chunks
			case 0x4110:
				if (!read_vertex_block(verts)) return 0;
				break;
				// TRI_MAPPINGCOORS: Vertices list
				// Chunk Length: 1 x unsigned short (# mapping points) + 2 x float (mapping coordinates) x (# mapping points) + sub chunks
			case 0x4140:
				if (!read_mapping_block(verts)) return 0;
				break;
			default:
				skip_chunk(chunk_len);
			} // end switch
		} // end while
		assert(ftell(fp) == end_pos);
		vector<counted_normal> normals;
		if (use_vertex_normals) {normals.resize(verts.size());}

		// build vertex lists and compute face normals
		// Note: we don't support vertices that are shared between faces of different materials
		for (vector<face_t>::const_iterator i = faces.begin(); i != faces.end(); ++i) {
			point pts[3];
			
			for (unsigned n = 0; n < 3; ++n) {
				unsigned const ix(i->ix[n]);
				assert(ix < verts.size());
				pts[n] = verts[n].v;
			}
			if (use_vertex_normals) {
				vector3d const normal(get_poly_norm(pts));
				UNROLL_3X(normals[i->ix[i_]].add_normal(normal);)
			}
		}
		model3d::proc_counted_normals(normals); // if use_vertex_normals

		// assign materials to faces
		for (face_mat_map_t::const_iterator i = face_materials.begin(); i != face_materials.end(); ++i) {
			for (vector<unsigned short>::const_iterator f = i->second.begin(); f != i->second.end(); ++f) {
				assert(*f < faces.size());
				assert(faces[*f].mat == -1); // material not yet assigned
				faces[*f].mat = i->first;
			}
		}
		vector<unsigned short> &def_mat(face_materials[-1]); // create default material

		for (unsigned i = 0; i < faces.size(); ++i) {
			if (faces[i].mat == -1) {def_mat.push_back(i);} // faces not assigned to a material get the default material
		}

		// add triangles to model for each material
		polygon_t tri;
		tri.resize(3);

		for (face_mat_map_t::const_iterator i = face_materials.begin(); i != face_materials.end(); ++i) {
			vntc_map_t vmap(0);

			for (vector<unsigned short>::const_iterator f = i->second.begin(); f != i->second.end(); ++f) {
				unsigned short *ixs(faces[*f].ix);
				point pts[3];
				UNROLL_3X(pts[i_] = verts[ixs[i_]].v;)
				vector3d const face_n(get_poly_norm(pts));

				for (unsigned j = 0; j < 3; ++j) {
					unsigned const ix(ixs[j]);
					vector3d const normal((!use_vertex_normals || (face_n != zero_vector && !normals[ix].is_valid())) ? face_n : normals[ix]);
					tri[j] = vert_norm_tc(pts[j], normal, verts[ix].t[0], verts[ix].t[1]);
				}
				model.add_triangle(tri, vmap, i->first, obj_id);
			}
		} // for i
		++obj_id;
		return 1;
	}

	bool read_and_proc_texture(unsigned chunk_len, int &tid, char const *const name) {
		string tex_fn;
		if (!read_texture(chunk_len, tex_fn)) {cerr << "Error reading texture " << name << endl; return 0;}
		check_and_bind(tid, tex_fn, 0, verbose);
		return 1;
	}

	bool read_material(unsigned read_len, material_t &cur_mat) {
		unsigned short chunk_id;
		unsigned chunk_len;
		long const end_pos(get_end_pos(read_len));

		while (ftell(fp) < end_pos) { // read each chunk from the file 
			if (!read_chunk_header(chunk_id, chunk_len)) return 0;

			switch (chunk_id) {
			case 0xA000: // material name
				if (!read_null_term_string(cur_mat.name)) return 0;
				break;
			case 0xA010: // material ambient color
				if (!read_color(cur_mat.ka)) return 0;
				break;
			case 0xA020: // material diffuse color
				if (!read_color(cur_mat.kd)) return 0;
				break;
			case 0xA030: // material specular color
				if (!read_color(cur_mat.ks)) return 0;
				break;
			case 0xA040: // material shininess
				if (!read_percentage(chunk_len, cur_mat.ns)) return 0;
				// FIXME: wrong (0.0003)?
				break;
			case 0xA050: // material transparency
				if (!read_percentage(chunk_len, cur_mat.alpha)) return 0;
				cur_mat.alpha = 1.0 - cur_mat.alpha; // convert from transparency to opacity
				break;
			case 0xA200: // texture map 1
				if (!read_and_proc_texture(chunk_len, cur_mat.d_tid, "texture map 1")) return 0;
				break;
			case 0xA230: // bump map
				if (!read_and_proc_texture(chunk_len, cur_mat.bump_tid, "bump map")) return 0;
				break;
			case 0xA220: // reflection map
				if (!read_and_proc_texture(chunk_len, cur_mat.refl_tid, "reflection map")) return 0;
				break;
			default:
				skip_chunk(chunk_len);
			} // end switch
		} // end while
		assert(ftell(fp) == end_pos);
		if (verbose) {cout << "Read material " << cur_mat.name << endl;}
		return 1;
	}

	bool read_texture(unsigned read_len, string &tex_name) {
		unsigned short chunk_id, map_tiling;
		unsigned chunk_len;
		long const end_pos(get_end_pos(read_len));

		while (ftell(fp) < end_pos) { // read each chunk from the file 
			if (!read_chunk_header(chunk_id, chunk_len)) return 0;

			switch (chunk_id) {
			case 0xA300: // mapping filename
				if (!read_null_term_string(tex_name)) return 0;
				break;
			case 0xA351: // mapping parameters
				if (!read_data(&map_tiling, sizeof(unsigned short), 1, "texture map tiling flags")) return 0;
				break;
			default:
				skip_chunk(chunk_len);
			} // end switch
		} // end while
		assert(ftell(fp) == end_pos);
		// FIXME: use or return map_tiling?
		return 1;
	}

public:
	file_reader_3ds_model(string const &fn, bool use_vertex_normals_, model3d &model_) :
	  file_reader_3ds(fn), model_from_file_t(fn, model_), use_vertex_normals(use_vertex_normals_), obj_id(0) {}

	bool read(geom_xform_t const &xf, bool verbose) {
		if (!file_reader_3ds::read(xf, verbose)) return 0;
		model.optimize(); // optimize vertices and remove excess capacity

		if (verbose) {
			cout << "bcube: "; model.get_bcube().print(); cout << endl;
			cout << "model stats: "; model.show_stats();
		}
		return 1;
	}
};


bool read_3ds_file_model(string const &filename, model3d &model, geom_xform_t const &xf, bool use_vertex_normals, bool verbose) {
	file_reader_3ds_model reader(filename, use_vertex_normals, model);
	return reader.read(xf, verbose);
}

bool read_3ds_file_pts(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, colorRGBA const &def_c, bool verbose) {
	file_reader_3ds_triangles reader(filename);
	return reader.read(ppts, xf, def_c, verbose);
}


