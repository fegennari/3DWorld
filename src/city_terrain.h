// 3D World - City Terrain and Road System
// by Frank Gennari
// 08/22/21

#pragma once

#include "city.h"
#include "mesh.h"

float const OUTSIDE_TERRAIN_HEIGHT  = 0.0;

struct rect_t {
	unsigned x1, y1, x2, y2;
	rect_t() : x1(0), y1(0), x2(0), y2(0) {}
	rect_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
	bool     is_valid() const {return (x1 < x2 && y1 < y2);}
	unsigned get_area() const {return (x2 - x1)*(y2 - y1);}
	bool operator== (rect_t const &r) const {return (x1 == r.x1 && y1 == r.y1 && x2 == r.x2 && y2 == r.y2);}
	bool has_overlap(rect_t const &r) const {return (x1 < r.x2 && y1 < r.y2 && r.x1 < x2 && r.y1 < y2);}
};

struct flatten_op_t : public rect_t {
	float z1, z2;
	bool dim=0;
	unsigned border=0, skip_six=0, skip_eix=0;
	flatten_op_t() : z1(0.0), z2(0.0) {}
	flatten_op_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_, float z1_, float z2_, bool dim_, unsigned border_) :
		rect_t(x1_, y1_, x2_, y2_), z1(z1_), z2(z2_), dim(dim_), border(border_) {}
};

class heightmap_query_t {
protected:
	float *heightmap;
	unsigned xsize, ysize;
public:
	flatten_op_t last_flatten_op;

	heightmap_query_t() : heightmap(nullptr), xsize(0), ysize(0) {}
	heightmap_query_t(float *hmap, unsigned xsize_, unsigned ysize_) : heightmap(hmap), xsize(xsize_), ysize(ysize_) {assert(heightmap != nullptr);}
	float get_x_value(int x) const {return get_xval(x - int(xsize)/2);} // convert from center to LLC
	float get_y_value(int y) const {return get_yval(y - int(ysize)/2);}
	int get_x_pos(float x) const {return (get_xpos(x) + int(xsize)/2);}
	int get_y_pos(float y) const {return (get_ypos(y) + int(ysize)/2);}
	float  get_height(unsigned x, unsigned y) const {return heightmap[y*xsize + x];} // Note: not bounds checked
	float &get_height(unsigned x, unsigned y)       {return heightmap[y*xsize + x];} // Note: not bounds checked
	float  get_height_bc(unsigned x, unsigned y) const {assert(x < xsize && y < ysize); return heightmap[y*xsize + x];}
	float  get_height_clamped(int x, int y) const {return (is_inside_terrain(x, y) ? get_height(x, y) : OUTSIDE_TERRAIN_HEIGHT);}
	bool is_normalized_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <  x2 && y1 <  y2 && x2 <= xsize && y2 <= ysize);}
	bool is_valid_region     (unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <= x2 && y1 <= y2 && x2 <= xsize && y2 <= ysize);}
	bool is_inside_terrain(int x, int y) const {return (x >= 0 && y >= 0 && x < (int)xsize && y < (int)ysize);}
	cube_t get_full_hmap_bcube() const {return get_cube_for_bounds(0, 0, xsize, ysize, 0.0);}
	cube_t get_cube_for_bounds(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) const;
	cube_t get_cube_for_cell(int x, int y) const;
	point  get_pos_for_cell (int x, int y) const {return point(get_x_value(x), get_y_value(y), get_height_clamped(x, y));}
	float get_height_at(float xval, float yval) const {return get_height_clamped(get_x_pos(xval), get_y_pos(yval));}
	float get_height_at(point const &pos      ) const {return get_height_at(pos.x, pos.y);}
	float get_road_zval_at_pt(point const &pos) const {return get_height_at(pos) + ROAD_HEIGHT;}
	bool any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2, bool check_border=0) const;
	void get_segment_end_pts(road_t const &r, unsigned six, unsigned eix, point &ps, point &pe) const;
	void flatten_region_to(cube_t const c, unsigned slope_width, bool decrease_only=0);
	void flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation, bool decrease_only=0);
	float flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, unsigned border,
		unsigned skip_six=0, unsigned skip_eix=0, bool stats_only=0, bool decrease_only=0, bridge_t *bridge=nullptr, tunnel_t *tunnel=nullptr);
	float flatten_for_road(road_t const &road, unsigned border, bool stats_only=0, bool decrease_only=0, bridge_t *bridge=nullptr, tunnel_t *tunnel=nullptr);
};

// structs/classes used for connector road path finding
struct road_endpoint_t {
	point pt;
	bool dim, dir;
	road_endpoint_t() : dim(0), dir(0) {}
	road_endpoint_t(point const &pt_, bool dim_, bool dir_) : pt(pt_), dim(dim_), dir(dir_) {}
};

struct road_cand_t {
	vector<point> pts;
	bool start_dim;
	float cost=0.0;
	road_cand_t(bool sdim=0) : start_dim(sdim) {}
	void clear() {pts.clear(); cost = 0.0;}
	bool valid() const {return !pts.empty();}
};

struct conn_isec_t {
	flatten_op_t fop;
	cube_t int_bcube;
	unsigned road_ix;
	bool dim, dir;
	conn_isec_t(heightmap_query_t const &hq, cube_t const &ibc, unsigned rix, bool dim_, bool dir_) : fop(hq.last_flatten_op), int_bcube(ibc), road_ix(rix), dim(dim_), dir(dir_) {}
};

struct city_road_connector_t {
	heightmap_query_t &hq; // Note: hq is not modifed
	vector<road_t> segments; // reused temporary
	rand_gen_t rgen;

	city_road_connector_t(heightmap_query_t &hq_) : hq(hq_) {}
	static bool get_closer_dir(cube_t const &A, cube_t const &B, bool dim);
	// roads
	float calc_road_cost(point const &p1, point const &p2);
	float calc_road_path_cost(vector<point> &pts);
	float find_route_between_points(point const &p1, point const &p2, vect_cube_t const &blockers, vector<point> &pts,
		cube_t const &bcube1, cube_t const &bcube2, float road_hwidth, bool dim1, bool dir1, bool dim2, bool dir2);
	bool segment_road(road_t const &road, bool check_only);
	// transmission lines
	bool is_tline_seg_valid(point const &p1, point const &p2, float max_ground_clearance) const;
	bool route_transmission_line(transmission_line_t &tline, vect_cube_t &blockers, float road_width, float road_spacing) const;
};

