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

	void read_vertex_block(geom_xform_t const &xf) {
		unsigned short num;
		fread(&num, sizeof(unsigned short), 1, fp);
		assert(verts.empty()); // can only get here once
		verts.resize(num);
		if (verbose) {printf("Number of vertices: %d\n", num);}
		for (int i = 0; i < num; i++) {
			fread(&verts[i].v.x, sizeof(float), 3, fp);
			xf.xform_pos(verts[i].v);
		}
	}

	void read_mapping_block() {
		unsigned short num;
		fread(&num, sizeof(unsigned short), 1, fp);
		assert(num == verts.size()); // must be one for each vertex
		if (verbose) {printf("Number of tcs: %d\n", num);}
		for (int i = 0; i < num; i++) {fread(verts[i].t, sizeof(float), 2, fp);}
	}

	virtual bool proc_other_chunks(unsigned chunk_id, unsigned chunk_len) {
		return 0;
	}

public:
	file_reader_3ds(string const &fn) : base_file_reader(fn) {}
	virtual ~file_reader_3ds() {}

	bool read(geom_xform_t const &xf, bool verbose_) {
		RESET_TIME;
		verbose = verbose_;
		if (!open_file(1)) return 0; // binary file
		cout << "Reading 3DS file " << filename << endl;
		unsigned short chunk_id;

		// FIXME: should really check all these fread return values
		while (fread(&chunk_id, 2, 1, fp) == 1) { // read each chunk from the file 
			//if (verbose) {printf("ChunkID: %x\n", chunk_id);}
			unsigned chunk_length;
			fread(&chunk_length, 4, 1, fp);
			//if (verbose) {printf("ChunkLength: %d\n", chunk_length);}

			switch (chunk_id) {
				// MAIN3DS: Main chunk, contains all the other chunks; length: 0 + sub chunks
			case 0x4d4d:
				break;
				// EDIT3DS: 3D Editor chunk, objects layout info; length: 0 + sub chunks
			case 0x3d3d:
				break;
				// EDIT_OBJECT: Object block, info for each object; length: len(object name) + sub chunks
			case 0x4000:
				for (int i = 0; ; ++i) {
					char l_char;
					fread(&l_char, 1, 1, fp);
					name.push_back(l_char); // add the null terminator or not?
					if (l_char == '\0') break;
				}
				break;
				// OBJ_TRIMESH: Triangular mesh, contains chunks for 3d mesh info; length: 0 + sub chunks
			case 0x4100:
				break;
				// TRI_VERTEXL: Vertices list
				// Chunk Length: 1 x unsigned short (# vertices) + 3 x float (vertex coordinates) x (# vertices) + sub chunks
			case 0x4110:
				read_vertex_block(xf);
				break;
				// TRI_MAPPINGCOORS: Vertices list
				// Chunk Length: 1 x unsigned short (# mapping points) + 2 x float (mapping coordinates) x (# mapping points) + sub chunks
			case 0x4140:
				read_mapping_block();
				break;

			default: // send to derived class reader, and skip chunk if it isn't handled
				if (!proc_other_chunks(chunk_id, chunk_length)) {fseek(fp, chunk_length-6, SEEK_CUR);}
			} // end switch
		}
		return 1;
	}
};


// ************************************************


class file_reader_3ds_triangles : public file_reader_3ds {

	colorRGBA def_color;
	vector<coll_tquad> *ppts;

	virtual bool proc_other_chunks(unsigned chunk_id, unsigned chunk_len) {
		assert(ppts != nullptr);
		unsigned short num;

		switch (chunk_id) {
			// TRI_FACEL1: Polygons (faces) list
			// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
		case 0x4120:
			fread(&num, sizeof(unsigned short), 1, fp);
			if (verbose) {printf("Number of polygons: %d\n", num);}

			for (int i = 0; i < num; i++) {
				unsigned short ix[3];
				fread(ix, sizeof(unsigned short), 3, fp);
				unsigned short l_face_flags; // ignored
				fread(&l_face_flags, sizeof(unsigned short), 1, fp);
				triangle tri;
				
				for (unsigned j = 0; j < 3; ++j) {
					assert(ix[j] < verts.size());
					tri.pts[j] = verts[ix[j]].v;
				}
				ppts->push_back(coll_tquad(tri, def_color));
			}
			return 1;
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

	model3d &model;
	vntc_map_t vmap;

	virtual bool proc_other_chunks(unsigned chunk_id, unsigned chunk_len) {
		unsigned short num;

		switch (chunk_id) {
			// TRI_FACEL1: Polygons (faces) list
			// Chunk Length: 1 x unsigned short (# polygons) + 3 x unsigned short (polygon points) x (# polygons) + sub chunks
		case 0x4120:
			fread(&num, sizeof(unsigned short), 1, fp);
			if (verbose) {printf("Number of polygons: %d\n", num);}
			polygon_t tri;
			tri.resize(3);
			int const mat_id(-1); // undefined
			unsigned const obj_id(0); // default
			vector<unsigned short> ixs(3*num);
			vector<vector3d> normals(verts.size());
			// FIXME: recalc_normals

			for (int i = 0; i < num; i++) {
				unsigned short *ix(&ixs.front() + 3*i);
				fread(ix, sizeof(unsigned short), 3, fp);
				unsigned short l_face_flags; // ignored
				fread(&l_face_flags, sizeof(unsigned short), 1, fp);
				point pts[3];

				for (unsigned j = 0; j < 3; ++j) {
					assert(ix[j] < verts.size());
					pts[j] = verts[ix[j]].v;
				}
				vector3d const normal(get_poly_norm(pts));
				UNROLL_3X(normals[ix[i_]] += normal;)
			}
			for (vector<vector3d>::iterator i = normals.begin(); i != normals.end(); ++i) {
				i->normalize();
			}
			for (int i = 0; i < num; i++) {
				for (unsigned j = 0; j < 3; ++j) {
					unsigned const ix(ixs[3*i + j]);
					tri[j] = vert_norm_tc(verts[ix].v, normals[ix], verts[ix].t[0], verts[ix].t[1]);
				}
				model.add_triangle(tri, vmap, mat_id, obj_id);
			}
			return 1;
		} // end switch
		return 0;
	}

public:
	file_reader_3ds_model(string const &fn, model3d &model_) : file_reader_3ds(fn), model(model_), vmap(0) {}

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


bool read_3ds_file_model(string const &filename, model3d &model, geom_xform_t const &xf, bool verbose) {
	file_reader_3ds_model reader(filename, model);
	return reader.read(xf, verbose);
}

bool read_3ds_file_pts(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, colorRGBA const &def_c, bool verbose) {
	file_reader_3ds_triangles reader(filename);
	return reader.read(ppts, xf, def_c, verbose);
}


