// 3D World - Navigation Grid class for use with building and city AI navigation
// by Frank Gennari
// 2-24-24
#pragma once

#include "pedestrians.h" // for ai_path_t

class cube_nav_grid {
protected:
	float radius=0.0;
	cube_t bcube, grid_bcube;
	unsigned num[2] = {};
	float   step[2] = {};
	bool invalid=0;
	uint8_t exclude_val=255; // set to a large value
	vector<uint8_t> nodes; // num[0] x num[1]; 0=open, >0=blocked
	vect_cube_t blockers_exp;

	struct ix_pair_t {
		unsigned x, y;
		ix_pair_t(unsigned x_, unsigned y_) : x(x_), y(y_) {}
		bool operator<(ix_pair_t const &p) const {return ((y == p.y) ? (x < p.x) : (y < p.y));} // needed for priority_queue
	};
	struct a_star_node_state_t {
		int16_t  came_from[2] = {-1,-1};
		uint16_t xy       [2] = { 0, 0};
		float g_score=0.0, f_score=0.0;
		void set(unsigned from_x, unsigned from_y, unsigned x, unsigned y) {came_from[0] = from_x; came_from[1] = from_y; xy[0] = x; xy[1] = y;}
	};
	static bool pt_contained_xy(point const &pt, vect_cube_t const &cubes);
	bool    are_ixs_valid(unsigned x, unsigned y) const {return (x < num[0] && y < num[1]);} // negative numbers will wrap around and still fail
	unsigned get_node_ix (unsigned x, unsigned y) const {assert(are_ixs_valid(x, y)); return (x + y* num[0]);}
	point    get_grid_pt (unsigned x, unsigned y, float zval) const {return point((grid_bcube.x1() + x*step[0]), (grid_bcube.y1() + y*step[1]), zval);}
	point    get_grid_pt (unsigned x, unsigned y) const {return get_grid_pt(x, y, grid_bcube.z1());}
	uint8_t  get_node_val(unsigned x, unsigned y) const {return nodes[get_node_ix(x, y)];}
	unsigned get_node_ix(point p) const;
	bool is_blocked(uint8_t val)            const {return (val > 0 && val != exclude_val);}
	bool is_blocked(unsigned x, unsigned y) const {return is_blocked(get_node_val(x, y));}
	void get_grid_ix_fp(point p, float gxy[2]) const;
	float get_distance(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const;
	bool find_open_node_closest_to(point const &p, point const &dest, unsigned &nx, unsigned &ny) const;
	virtual bool check_line_intersect(point const &p1, point const &p2, float radius) const;
	void get_region_xy_bounds(cube_t const &region, unsigned &x1, unsigned &x2, unsigned &y1, unsigned &y2) const;
	void make_region_walkable(cube_t const &region);
public:
	virtual ~cube_nav_grid() {}
	bool is_built() const {return !bcube.is_all_zeros();} // can't test on !nodes.empty() in case the room is too small to have any nodes
	bool is_valid() const {return (!invalid && is_built());}
	void invalidate() {invalid = 1;}
	void build(cube_t const &bcube_, vect_cube_t const &blockers, float radius_, bool add_edge_pad, bool no_blocker_expand);
	bool find_path(point const &p1, point const &p2, ai_path_t &path) const;
};
