// 3D World - 3D Studio Max File Reader
// by Frank Gennari
// 11/1/2014
// Reference: http://www.spacesimulator.net/wiki/index.php?title=Tutorials:3ds_Loader

#include "3DWorld.h"
#include "model3d.h"
#include "file_reader.h"


class file_reader_3ds : public base_file_reader {

protected:
	string name; // unused?
	vector<vert_tc_t> verts;

	bool read_data(void *dest, size_t sz, size_t count, char const *const str) {
		size_t const nread(fread(dest, sz, count, fp));
		if (nread == count) return 1;
		cerr << "Error reading 3DS file " << str << " data: expected " << count << " elements of size " << sz << " but got " << nread << endl;
		return 0;
	}

	bool read_chunk_header(unsigned short &chunk_id, unsigned &chunk_len) {
		if (!read_data(&chunk_id, 2, 1, "chunk id")) return 0;
		if (verbose) {printf("ChunkID: %x\n", chunk_id);}
		if (!read_data(&chunk_len, 4, 1, "chunk length")) return 0;
		if (verbose) {printf("ChunkLength: %d\n", chunk_len);}
		return 1;
	}

	bool read_vertex_block(geom_xform_t const &xf) {
		unsigned short num;
		if (!read_data(&num, sizeof(unsigned short), 1, "vertex size")) return 0;
		//assert(verts.empty()); // can only get here once?
		verts.resize(num);
		if (verbose) {printf("Number of vertices: %d\n", num);}

		for (int i = 0; i < num; i++) {
			if (!read_data(&verts[i].v.x, sizeof(float), 3, "vertex xyz")) return 0;
			xf.xform_pos(verts[i].v);
		}
		return 1;
	}

	bool read_mapping_block() {
		unsigned short num;
		if (!read_data(&num, sizeof(unsigned short), 1, "mapping size")) return 0;
		assert(num == verts.size()); // must be one for each vertex
		if (verbose) {printf("Number of tex coords: %d\n", num);}

		for (int i = 0; i < num; i++) {
			if (!read_data(verts[i].t, sizeof(float), 2, "mapping uv")) return 0;
		}
		return 1;
	}

	bool read_one_face(unsigned short ix[3]) {
		if (!read_data(ix, sizeof(unsigned short), 3, "face indices")) return 0;
		UNROLL_3X(assert(ix[i_] < verts.size());)
		unsigned short l_face_flags; // 3 LSB are edge orderings (AB BC CA): ignored
		if (!read_data(&l_face_flags, sizeof(unsigned short), 1, "face flags")) return 0;
		//cout << "flags: " << l_face_flags << endl;
		return 1;
	}

	bool read_null_term_string(string &str) {
		while (1) {
			char l_char;
			if (!read_data(&l_char, 1, 1, "string char")) return 0;
			str.push_back(l_char); // add the null terminator or not?
			if (l_char == '\0') return 1;
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

	virtual bool proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {return 0;}

public:
	file_reader_3ds(string const &fn) : base_file_reader(fn) {}
	virtual ~file_reader_3ds() {}

	bool read(geom_xform_t const &xf, bool verbose_) {
		RESET_TIME;
		verbose = verbose_;
		if (!open_file(1)) return 0; // binary file
		cout << "Reading 3DS file " << filename << endl;
		unsigned short chunk_id;
		unsigned chunk_len;

		while (read_chunk_header(chunk_id, chunk_len)) { // read each chunk from the file 
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
				// OBJ_TRIMESH: Triangular mesh, contains chunks for 3d mesh info; length: 0 + sub chunks
			case 0x4100:
				break;
				// TRI_VERTEXL: Vertices list
				// Chunk Length: 1 x unsigned short (# vertices) + 3 x float (vertex coordinates) x (# vertices) + sub chunks
			case 0x4110:
				if (!read_vertex_block(xf)) return 0;
				break;
				// TRI_MAPPINGCOORS: Vertices list
				// Chunk Length: 1 x unsigned short (# mapping points) + 2 x float (mapping coordinates) x (# mapping points) + sub chunks
			case 0x4140:
				if (!read_mapping_block()) return 0;
				break;

			default: // send to derived class reader, and skip chunk if it isn't handled
				if (!proc_other_chunks(chunk_id, chunk_len)) {fseek(fp, chunk_len-6, SEEK_CUR);}
			} // end switch
		}
		return 1;
	}
};


// ************************************************


class file_reader_3ds_triangles : public file_reader_3ds {

	colorRGBA def_color;
	vector<coll_tquad> *ppts;

	virtual bool proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {
		assert(ppts != nullptr);
		unsigned short num;

		switch (chunk_id) {
			// TRI_FACEL1: Polygons (faces) list
			// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
		case 0x4120:
			if (!read_data(&num, sizeof(unsigned short), 1, "number of faces")) return 0;
			if (verbose) {printf("Number of polygons: %d\n", num);}

			for (int i = 0; i < num; i++) {
				unsigned short ix[3];
				read_one_face(ix);
				triangle tri;
				UNROLL_3X(tri.pts[i_] = verts[ix[i_]].v;)
				ppts->push_back(coll_tquad(tri, def_color));
			}
			return 1; // handled
		} // end switch
		return 0;
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


class file_reader_3ds_model : public file_reader_3ds {

	bool use_vertex_normals;
	model3d &model;

	virtual bool proc_other_chunks(unsigned short chunk_id, unsigned chunk_len) {

		switch (chunk_id) {
			// TRI_FACEL1: Polygons (faces) list
			// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
		case 0x4120:
		{
			unsigned short num;
			if (!read_data(&num, sizeof(unsigned short), 1, "number of faces")) return 0;
			if (verbose) {printf("Number of polygons: %d\n", num);}
			polygon_t tri;
			tri.resize(3);
			int const mat_id(-1); // undefined
			unsigned const obj_id(0); // default
			vector<unsigned short> ixs(3*num);
			vector<counted_normal> normals;
			if (use_vertex_normals) {normals.resize(verts.size());}

			// read faces, build vertex lists, and compute face normals
			for (int i = 0; i < num; i++) {
				unsigned short *ix(&ixs.front() + 3*i);
				read_one_face(ix);
				swap(ix[0], ix[2]); // reverse triangle vertex ordering to agree with 3DWorld coordinate system
				point pts[3];
				UNROLL_3X(pts[i_] = verts[ix[i_]].v;)

				if (use_vertex_normals) {
					vector3d const normal(get_poly_norm(pts));
					UNROLL_3X(normals[ix[i_]].add_normal(normal);)
				}
			}
			model3d::proc_counted_normals(normals); // if use_vertex_normals

			// compute vertex normals and add triangles to the model
			vntc_map_t vmap(0);

			for (int i = 0; i < num; i++) {
				point pts[3];
				UNROLL_3X(pts[i_] = verts[ixs[3*i + i_]].v;)
				vector3d const face_n(get_poly_norm(pts));

				for (unsigned j = 0; j < 3; ++j) {
					unsigned const ix(ixs[3*i + j]);
					vector3d const normal((!use_vertex_normals || (face_n != zero_vector && !normals[ix].is_valid())) ? face_n : normals[ix]);
					tri[j] = vert_norm_tc(pts[j], normal, verts[ix].t[0], verts[ix].t[1]);
				}
				model.add_triangle(tri, vmap, mat_id, obj_id);
			}
			return 1; // handled
		}

		// faces data
		case 0x4130: // Faces Material
			break;
		case 0x4150: // Smoothing Group List
			// nfaces*4bytes: Long int where the nth bit indicates if the face belongs to the nth smoothing group
			break;

		case 0xAFFF: // material
			{
				return 0; // FIXME - incomplete
				//material_t cur_mat; // FIXME
				//return read_material(chunk_len, cur_mat); // handled
			}
		} // end switch
		return 0;
	}

	bool proc_texture(string const &filename, int &tid) {
		// FIXME
		return 0;
	}

	bool read_material(unsigned read_len, material_t &cur_mat) {
		unsigned short chunk_id;
		unsigned chunk_len;
		string tex_fn;

		// FIXME: exit loop using read_len and ftell()
		while (read_chunk_header(chunk_id, chunk_len)) { // read each chunk from the file 
			switch (chunk_id) {
			case 0xA000: // material name
				if (!read_null_term_string(cur_mat.name)) return 0;
			case 0xA010: // material ambient color
				if (!read_color(cur_mat.ka)) return 0;
				break;
			case 0xA020: // material diffuse color
				if (!read_color(cur_mat.kd)) return 0;
				break;
			case 0xA030: // material specular color
				if (!read_color(cur_mat.ks)) return 0;
				break;
			case 0xA200: // texture map 1
				if (!read_texture(chunk_len, tex_fn)) return 0;
				if (!proc_texture(tex_fn, cur_mat.d_tid)) return 0;
				break;
			case 0xA230: // bump map
				if (!read_texture(chunk_len, tex_fn)) return 0;
				if (!proc_texture(tex_fn, cur_mat.bump_tid)) return 0;
				break;
			case 0xA220: // reflection map
				if (!read_texture(chunk_len, tex_fn)) return 0;
				if (!proc_texture(tex_fn, cur_mat.refl_tid)) return 0;
				break;
			default: return 0;
			} // end switch
		} // end while
		return 1;
	}

	bool read_texture(unsigned read_len, string &tex_name) {
		unsigned short chunk_id;
		unsigned chunk_len;

		// FIXME: exit loop using read_len and ftell()
		while (read_chunk_header(chunk_id, chunk_len)) { // read each chunk from the file 
			switch (chunk_id) {
			case 0xA300: // mapping filename
				if (!read_null_term_string(tex_name)) return 0;
				break;
			case 0xA351: // mapping parameters
				// FIXME: WRITE
				break;
			default: return 0;
			} // end switch
		} // end while
		return 1;
	}

public:
	file_reader_3ds_model(string const &fn, bool use_vertex_normals_, model3d &model_) :
	  file_reader_3ds(fn), use_vertex_normals(use_vertex_normals_), model(model_) {}

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


