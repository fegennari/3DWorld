// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/29/05

#include "3DWorld.h"
#include "shape_line3d.h"
#include "collision_detect.h"
#include "draw_utils.h"
#include <fstream>


extern vector3d up_norm;


// ************** SHAPE3D *************


bool shape3d::alloc_shape(unsigned npoints, unsigned nfaces, unsigned ncolors) {

	if (npoints < 3 || nfaces < 1) return 0;
	points.resize(npoints);
	faces.resize(nfaces);
	colors.resize(ncolors);
	return 1;
}


// similar to the object file reader
bool shape3d::read_from_file(char *filename) {

	float x, y, z;
	char  letter;
	unsigned ix, iy, iz, nfaces(0), nverts(0), ncolors(0), color_id;
	colorRGBA color;
	assert(filename);
	FILE *fp;
	if (!open_file(fp, filename, "shape3d")) return 0;

	if (fscanf(fp, "%c", &letter) == 1 && letter == 'n') {	
		if (fscanf(fp, "%i%i%i\n", &ncolors, &nverts, &nfaces) != 3) {
			cout << "Error reading number of vertices and faces from file '" << filename << "'." << endl;
			fclose(fp);
			return 0;
		}
	}
	else { // Count the number of vertices, faces, and colors
		ungetc((int)letter, fp);
		
		while(!feof(fp)) {
			if (fscanf(fp, "%c%f%f%f", &letter, &x, &y, &z) != 4) {
				cout << "Error reading entry from file '" << filename << "'." << endl;
				fclose(fp);
				return 0;
			}
			switch (letter) {
			case 'v': nverts++;  break;
			case 'f': nfaces++;  break;
			case 'c': ncolors++; break;
			default:
				cout << "Error: Invalid letter in input file: " << letter << endl;
				fclose(fp);
				return 0;
			}
		}
		fclose(fp);
		fp = fopen(filename, "r");
	}
	assert(fp);
	
	if (nverts < 3 || nfaces < 1) {
		cout << "Error: Mesh in file '" << filename << "' must have at least three vertices and one face: " << nverts << ", " << nfaces << "." << endl;
		fclose(fp);
		return 0;
	}
	assert(alloc_shape(nverts, nfaces, ncolors));

	// Read the colors
	for(unsigned i = 0; i < ncolors; i++) {
		if (fscanf(fp, "%c%f%f%f%f%f%f%i\n", &letter, &color.R, &color.G, &color.B, &color.A, &colors[i].spec1, &colors[i].spec2, &colors[i].tid) != 8) {
			cout << "Error reading color from file '" << filename << "'." << endl;
			fclose(fp);
			return 0;
		}
		if (letter != 'c') {
			cout << "Error: Expecting color " << i << " in input file but got character " << letter << " <" << int((unsigned char)letter) << ">" << endl;
			fclose(fp);
			return 0;
		}
		colors[i].c = color;
	}

	// Read the vertices
	for(unsigned i = 0; i < nverts; i++) {
		if (fscanf(fp, "%c%f%f%f\n", &letter, &points[i].x, &points[i].y, &points[i].z) != 4) {
			cout << "Error reading vertex from file '" << filename << "'." << endl;
			fclose(fp);
			return 0;
		}
		if (letter != 'v') { // in case vertices and faces are mixed together
			cout << "Error: Expecting vertex " << i << " in input file but got character " << letter << " <" << int((unsigned char)letter) << ">" << endl;
			fclose(fp);
			return 0;
		}
	}

	// Read the faces
	for(unsigned i = 0; i < nfaces; i++) {
		if (fscanf(fp, "%c%u%u%u%u\n", &letter, &ix, &iy, &iz, &color_id) != 5) {
			cout << "Error reading face from file '" << filename << "'." << endl;
			fclose(fp);
			return 0;
		}
		if (letter != 'f') { // in case vertices and faces are mixed together
			cout << "Error: Expecting face " << i << " in input file." << endl;
			fclose(fp);
			return 0;
		}
		if (color_id >= max(1U, ncolors)) {
			cout << "Illegal type: " << color_id << "." << endl;
			fclose(fp);
			return 0;
		}
		if (ix < 1 || iy < 1 || iz < 1 || ix > nverts || iy > nverts || iz > nverts) {
			cout << "Error: Face " << i << " references nonexistant vertex (" << ix << ", " << iy << ", " << iz << ")." << endl;
			fclose(fp);
			return 0;
		}
		faces[i].v[2]     = ix-1;
		faces[i].v[1]     = iy-1;
		faces[i].v[0]     = iz-1;
		faces[i].color_id = color_id;
	}
	fclose(fp);
	return 1;
}


void shape3d::gen_face_normals() {

	for (unsigned i = 0; i < faces.size(); ++i) {
		get_face_normal(i);
	}
}


void shape3d::get_face_normal(unsigned face_id) {

	assert(face_id < faces.size());
	unsigned const *const verts(faces[face_id].v);
	assert(verts != NULL);
	assert(verts[0] < points.size() && verts[1] < points.size() && verts[2] < points.size());
	get_normal(points[verts[0]], points[verts[1]], points[verts[2]], faces[face_id].norm, 1);
}


void shape3d::get_triangle_center(point &center, unsigned face_id, unsigned quality) {

	assert(face_id < faces.size());
	unsigned const *const verts(faces[face_id].v);
	assert(verts != NULL);
	assert(verts[0] < points.size() && verts[1] < points.size() && verts[2] < points.size());
	point const p1(points[verts[0]]), p2(points[verts[1]]), p3(points[verts[2]]);

	// faster than above but not completely accurate
	if (quality == 0) {
		for (unsigned i = 0; i < 3; ++i) {
			center[i] = (p1[i] + p2[i] + p3[i])/3.0;
		}
		return;
	}
	point m1, m2;

	for (unsigned i = 0; i < 3; ++i) {
		m1[i] = 0.5*(p1[i] + p2[i]);
		m2[i] = 0.5*(p1[i] + p3[i]);
	}
	vector3d const v1(vector3d(m1, p3).get_norm()), v2(vector3d(m2, p2).get_norm());
	vector3d const v12(cross_product(v1, v2).get_norm()), vp2p3(p2, p3);
	float const s1(vector_determinant(vp2p3, v2, v12));///magsq12;
	float const s2(vector_determinant(vp2p3, v1, v12));///magsq12;

	for (unsigned i = 0; i < 3; ++i) {
		center[i] = 0.5*(p3[i] + v1[i]*s1 + p2[i] + v2[i]*s2);
	}
}


void shape3d::add_vertex(unsigned vertex, unsigned face_id, unsigned &face_counter) {

	assert(face_id < faces.size());
	unsigned *verts(faces[face_id].v);
	assert(verts != NULL);
	unsigned const p0(verts[0]), p1(verts[1]), p2(verts[2]);
	assert(p0 < points.size() && p1 < points.size() && p2 < points.size());
	verts[2] = vertex;
	verts    = faces[face_counter++].v;
	verts[0] = p2;
	verts[1] = p0;
	verts[2] = vertex;
	verts    = faces[face_counter++].v;
	verts[0] = p1;
	verts[1] = p2;
	verts[2] = vertex;
}


void shape3d::draw(bool skip_color_set) const {

	if (points.empty()) return;
	
	if (points.size() < 3 || faces.empty()) {
		cout << "Invalid shape: Cannot draw." << endl;
		return;
	}
	if (!skip_color_set) {
		enable_blend();
		set_color(color);
	}
	unsigned lcid(0);
	if (colors.empty() && tid >= 0) select_texture(tid);
	vector<vert_norm_tc> verts;

	for (unsigned i = 0; i < faces.size(); ++i) {
		if (!colors.empty()) {
			unsigned const color_id(faces[i].color_id);

			if (i == 0 || color_id != lcid) {
				draw_and_clear_verts(verts, GL_TRIANGLES);
				select_texture(colors[color_id].tid);
				if (!skip_color_set) {set_color(colors[color_id].c);}
				set_specular(colors[color_id].spec1, colors[color_id].spec2);
				lcid = color_id;
			}
		}
		int const max_dim(get_max_dim(faces[i].norm));

		for (unsigned j = 0; j < 3; ++j) {
			unsigned const index(faces[i].v[j]);
			assert(index < points.size());
			point const p(points[index]);
			int const d1[3] = {1,0,0}, d2[3] = {2,2,1};
			verts.push_back(vert_norm_tc((p*scale + pos), faces[i].norm, tex_scale*p[d1[max_dim]], tex_scale*p[d2[max_dim]]));
		}
	}
	draw_verts(verts, GL_TRIANGLES);
	if (!skip_color_set) {disable_blend();}
	if (tid >= 0) {glDisable(GL_TEXTURE_2D);}
	if (!colors.empty()) {set_specular(0.0, 0.0);}
}


void shape3d::get_triangle_verts(vector<vert_norm_tc> &verts) const {

	if (points.empty()) return;
	assert(points.size() >= 3 && !faces.empty());
	assert(colors.empty()); // not supported

	for (unsigned i = 0; i < faces.size(); ++i) {
		int const max_dim(get_max_dim(faces[i].norm));

		for (unsigned j = 0; j < 3; ++j) {
			unsigned const index(faces[i].v[j]);
			assert(index < points.size());
			point const p(points[index]);
			int const d1[3] = {1,0,0}, d2[3] = {2,2,1};
			verts.push_back(vert_norm_tc((p*scale + pos), faces[i].norm, tex_scale*p[d1[max_dim]], tex_scale*p[d2[max_dim]]));
		}
	}
}


void shape3d::add_cobjs(vector<int> &cids, bool draw) {

	point points2[3];

	for (unsigned i = 0; i < faces.size(); ++i) {
		for (unsigned j = 0; j < 3; ++j) {
			points2[j] = points[faces[i].v[j]]*scale + pos;
		}
		scolor const &sc(colors[faces[i].color_id]);
		cids.push_back(add_coll_polygon(points2, 3, cobj_params(0.7, sc.c, draw, 0, NULL, 0, sc.tid), 0.0));
	}
}


void shape3d::destroy() {

	points.clear();
	faces.clear();
	colors.clear();
	color = BLACK; // reset in case it is to be used again
	pos   = all_zeros;
	scenery_obj::destroy();
}


// ************** LINE3D *************


void line3d::draw(bool draw_as_tquads) const {

	if (points.empty()) return;
	
	if (points.size() < 2) {
		cout << "Invalid line: Cannot draw." << endl;
		return;
	}
	if (draw_as_tquads) {
		float const w(0.01*width);
		line_tquad_draw_t drawer;

		for (unsigned i = 1; i < points.size(); ++i) {
			drawer.add_line_tquad(points[i-1], points[i], w, w, color, color,
				((i > 1) ? &points[i-2] : NULL), ((i+1 < points.size()) ? &points[i+1] : NULL));
		}
		drawer.draw(GL_QUADS);
	}
	else {
		//glEnable(GL_BLEND);
		//glEnable(GL_LINE_SMOOTH);
		glLineWidth(width);
		glDisable(GL_LIGHTING);
		color.do_glColor();
		vector<vert_wrap_t> verts;
		verts.reserve(2*(points.size() - 1));

		for (unsigned i = 1; i < points.size(); ++i) {
			for (unsigned d = 0; d < 2; ++d) {verts.push_back(points[i-!d]);}
		}
		draw_verts(verts, GL_LINES);
		glLineWidth(1.0);
		glEnable(GL_LIGHTING);
		//glDisable(GL_BLEND);
		//glDisable(GL_LINE_SMOOTH);
	}
}


void line3d::destroy() {

	points.clear();
	color = BLACK;
}


