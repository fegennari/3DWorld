// 3D World - 3D Studio Max File Reader
// by Frank Gennari
// 11/1/2014
// Reference: http://www.spacesimulator.net/wiki/index.php?title=Tutorials:3ds_Loader

#include "3DWorld.h"
#include "model3d.h"
#include "file_reader.h"


class file_reader_3ds : public base_file_reader {

public:
	file_reader_3ds(string const &fn) : base_file_reader(fn) {}

	bool read(vector<coll_tquad> *ppts, geom_xform_t const &xf, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		cout << "Reading 3DS file " << filename << endl;
		// FIXME: WRITE
		return 1;
	}
};


// ************************************************


class file_reader_3ds_model : public file_reader_3ds {

	model3d &model;

public:
	file_reader_3ds_model(string const &fn, model3d &model_) : file_reader_3ds(fn), model(model_) {}

	bool read(geom_xform_t const &xf, bool recalc_normals, bool verbose) {
		RESET_TIME;
		if (!open_file()) return 0;
		cout << "Reading 3DS file " << filename << endl;
		// FIXME: WRITE
		return 1;
	}
};


bool read_3ds_file_model(string const &filename, model3d &model, geom_xform_t const &xf, bool recalc_normals, bool verbose) {
	file_reader_3ds_model reader(filename, model);
	return reader.read(xf, recalc_normals, verbose);
}

bool read_3ds_file_pts(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, bool verbose) {
	file_reader_3ds reader(filename);
	return reader.read(ppts, xf, verbose);
}


