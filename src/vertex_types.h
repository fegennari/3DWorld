// 3D World - Vertex Types
// by Frank Gennari
// 1/5/20

struct vert_norm { // size = 24
	point v;
	vector3d n;
	typedef vector3d normal_type;
	vert_norm() {}
	vert_norm(point const &v_, vector3d const &n_) : v(v_), n(n_) {}
	void assign(point const &v_, vector3d const &n_) {v = v_; n = n_;}
	void set_norm(vector3d const &n_) {n = n_;}
	vector3d const &get_norm() const {return n;}
	void invert_normal() {n = -n;}
	bool operator< (vert_norm const &p) const {return ((v == p.v) ? (n < p.n) : (v < p.v));}
	bool operator==(vert_norm const &p) const {return (v == p.v && n == p.n);}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void set_vbo_arrays_shadow(bool include_tcs);
	static void unset_attrs() {}
};


struct norm_comp { // size = 4
	char n[3];
	char pad; // unused padding
	norm_comp() : pad(0) {set_norm_to_zero();}
	norm_comp(vector3d const &n_) : pad(0) {set_norm(n_);}
	void set_norm_to_zero() {n[0] = n[1] = n[2] = 0;}
	void set_ortho_norm(unsigned dim, bool dir) {assert(dim < 3); set_norm_to_zero(); n[dim] = (dir ? 127 : -128);}
	void set_norm(norm_comp const &n_) {UNROLL_3X(n[i_] = n_.n[i_];)}
	void set_norm(vector3d const &n_) {UNROLL_3X(n[i_] = char(max(-128, min(127, int(127.0*n_[i_]))));)}
	void set_norm(char const n_[3]) {UNROLL_3X(n[i_] = n_[i_];)}
	void set_norm_no_clamp(vector3d const &n_) {UNROLL_3X(n[i_] = int(127.0*n_[i_]);)}
	vector3d get_norm() const {return vector3d(n[0]/127.0, n[1]/127.0, n[2]/127.0);}
	void invert_normal_dim(unsigned d) {assert(d < 3); n[d] = ((n[d] == -128) ? 127 : ((n[d] == 127) ? -128 : -n[d]));} // careful to not wraparound for -128 => 128
	void invert_normal() {UNROLL_3X(invert_normal_dim(i_);)}
};


// unused
struct norm_xy { // size = 8
	float x, y; // z can be calculated in the shader as sqrt(1 - x*x - y*y), as long as z >= 0.0
	norm_xy() : x(0.0f), y(0.0f) {}
	norm_xy(vector3d const &n) {set_norm(n);}
	void set_norm(vector3d const &n) {assert(n.z >= 0.0); x = n.x; y = n.y;} 
	void ensure_normalized_and_set(vector3d const &n) {assert(n.z >= 0.0); float const mag(n.mag()); x = n.x/mag; y = n.y/mag;}
	vector3d get_norm() const {return vector3d(x, y, sqrt(max(0.0f, (1.0f - x*x - y*y))));}
};


struct vert_wrap_t { // size = 12; so we can put the vertex first
	point v;
	vert_wrap_t() {}
	vert_wrap_t(point const &v_) : v(v_) {}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void unset_attrs() {}
};


struct vert_tc_t : public vert_wrap_t { // size = 20
	float t[2];
	typedef vert_tc_t non_color_class;
	vert_tc_t() {t[0] = t[1] = 0.0f;}
	vert_tc_t(point const &v_, float ts, float tt) : vert_wrap_t(v_) {t[0] = ts; t[1] = tt;}
	vert_tc_t(float x, float y, float z, float ts, float tt) : vert_wrap_t(point(x, y, z)) {t[0] = ts; t[1] = tt;}
	void assign(point const &v_, float ts, float tt) {v = v_; t[0] = ts; t[1] = tt;}
	bool operator==(vert_tc_t const &p) const {return (v == p.v && t[0] == p.t[0] && t[1] == p.t[1]);}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_comp : public vert_wrap_t, public norm_comp { // size = 16
	typedef norm_comp normal_type;
	vert_norm_comp() {}
	vert_norm_comp(vert_norm const &vn) : vert_wrap_t(vn.v), norm_comp(vn.n) {}
	vert_norm_comp(point const &v_, vector3d  const &n_) : vert_wrap_t(v_), norm_comp(n_) {}
	vert_norm_comp(point const &v_, norm_comp const &n_) : vert_wrap_t(v_), norm_comp(n_) {}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	void swap_dims(unsigned d1, unsigned d2) {assert(d1 < 3 && d2 < 3); swap(v[d1], v[d2]); swap(n[d1], n[d2]);}
	void invert_dim(unsigned d) {assert(d < 3); v[d] = -v[d]; invert_normal_dim(d);}
};


struct vert_norm_comp_tc : public vert_norm_comp { // size = 24
	float t[2];
	typedef vert_norm_comp_tc non_color_class;
	vert_norm_comp_tc() {t[0] = t[1] = 0.0f;}
	vert_norm_comp_tc(vert_norm_comp const &vn, float ts, float tt) : vert_norm_comp(vn) {t[0] = ts; t[1] = tt;}
	vert_norm_comp_tc(point const &v_, norm_comp const &n_, float ts, float tt) : vert_norm_comp(v_, n_) {t[0] = ts; t[1] = tt;}
	vert_norm_comp_tc(point const &v_, vector3d  const &n_, float ts, float tt) : vert_norm_comp(v_, n_) {t[0] = ts; t[1] = tt;}
	vert_norm_comp_tc(point const &v_, vector3d  const &n_, float const tc[2] ) : vert_norm_comp(v_, n_) {t[0] = tc[0]; t[1] = tc[1];}
	void set_tc(float tx, float ty) {t[0] = tx; t[1] = ty;}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_comp_tc_comp : public vert_norm_comp { // size = 20
	short t[2]; // could even use char
	vert_norm_comp_tc_comp() {t[0] = t[1] = 0;}
	vert_norm_comp_tc_comp(point const &v_, vector3d const &n_, float ts, float tt) : vert_norm_comp(v_, n_) {t[0] = 32767*ts; t[1] = 32767*tt;}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_tc : public vert_norm { // size = 32
	float t[2];
	typedef vert_norm_tc non_color_class;
	vert_norm_tc() {t[0] = t[1] = 0.0f;}
	vert_norm_tc(point const &v_, vector3d const &n_, float ts, float tt) : vert_norm(v_, n_) {t[0] = ts;    t[1] = tt;   }
	vert_norm_tc(point const &v_, vector3d const &n_, float const t_[2])  : vert_norm(v_, n_) {t[0] = t_[0]; t[1] = t_[1];}
	vert_norm_tc(vert_norm const &vn, float ts=0.0, float tt=1.0) : vert_norm(vn) {t[0] = ts; t[1] = tt;}
	void assign(point const &v_, vector3d const &n_, float ts, float tt) {v = v_; n = n_; t[0] = ts; t[1] = tt;}

	bool operator<(vert_norm_tc const &p) const {
		if (v < p.v) return 1;
		if (p.v < v) return 0;
		if (n < p.n) return 1;
		if (p.n < n) return 0;
		if (t[0] < p.t[0]) return 1;
		if (p.t[0] < t[0]) return 0;
		return (t[1] < p.t[1]);
	}
	bool operator==(vert_norm_tc const &p) const {return (v == p.v && n == p.n && t[0] == p.t[0] && t[1] == p.t[1]);}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void set_vbo_arrays_shadow(bool include_tcs);
};


struct vert_norm_tc_tan : public vert_norm_tc { // size = 48
	vector4d tangent;

	vert_norm_tc_tan() {}
	vert_norm_tc_tan(vert_norm_tc const &vntc, vector4d const &tangent_=vector4d(0,0,0,0))
		: vert_norm_tc(vntc), tangent(tangent_) {}
	vert_norm_tc_tan(point const &v_, vector3d const &n_, float ts, float tt, vector4d const &tangent_=vector4d(0,0,0,0))
		: vert_norm_tc(v_, n_, ts, tt), tangent(tangent_) {}

	bool operator<(vert_norm_tc_tan const &p) const {
		if (v < p.v) return 1;
		if (p.v < v) return 0;
		if (n < p.n) return 1;
		if (p.n < n) return 0;
		if (t[0] < p.t[0]) return 1;
		if (p.t[0] < t[0]) return 0;
		if (t[1] < p.t[1]) return 1;
		if (p.t[1] < t[1]) return 0;
		return (tangent < p.tangent);
	}
	//bool operator==(vert_norm_tc_tan const &p) const {}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void set_vbo_arrays_shadow(bool include_tcs);
	static void unset_attrs();
};


struct color_wrapper { // size = 4, can be used in a union
	unsigned char c[4]; // Note: c[3] (alpha component) is not used in all cases

	color_wrapper() {c[0] = c[1] = c[2] = c[3] = 0;}
	color_wrapper(colorRGBA const &c_) {set_c4(c_);}
	color_wrapper(colorRGB  const &c_) {set_c3(c_);}
	bool operator==(color_wrapper const &w) const {return (c[0] == w.c[0] && c[1] == w.c[1] && c[2] == w.c[2] && c[3] == w.c[3]);}
	template<typename T> void set_c3(T const &c_) {UNROLL_3X(c[i_] = (unsigned char)(255.0*CLIP_TO_01(c_[i_]));) c[3] = 255;}
	void set_c4(colorRGBA const &c_) {UNROLL_4X(c[i_]  = (unsigned char)(255.0*CLIP_TO_01(c_[i_]));)}
	void add_c4(colorRGBA const &c_) {UNROLL_4X(c[i_] += (unsigned char)(255.0*CLIP_TO_01(c_[i_]));)}
	void copy_color(color_wrapper const &cw) {UNROLL_4X(c[i_] = cw.c[i_];)}
	void copy_color(unsigned char const *const c_, bool has_alpha=0) {UNROLL_3X(c[i_] = c_[i_];) c[3] = (has_alpha ? c_[3] : 255);}
	colorRGB  get_c3() const {return colorRGB(c[0]/255.0, c[1]/255.0, c[2]/255.0);}
	colorRGBA get_c4() const {return colorRGBA(get_c3(), c[3]/255.0);}
	static bool is_compressed() {return 1;}
};

struct color_wrapper_ctor : public color_wrapper { // size = 4
	color_wrapper_ctor() {}
	color_wrapper_ctor(colorRGB  const &color) {set_c3(color);}
	color_wrapper_ctor(colorRGBA const &color) {set_c4(color);}
};


struct color_wrapper_float { // size = 16
	colorRGBA c; // Note: c[3] (alpha component) is not used in all cases

	template<typename T> void set_c3(T const &c_) {c = c_; c.A = 1.0;}
	void set_c4(colorRGBA const &c_) {c = c_;}
	colorRGB  get_c3() const {return colorRGB(c.R, c.G, c.B);}
	colorRGBA get_c4() const {return c;}
	static bool is_compressed() {return 0;}
};


struct vert_color : public color_wrapper { // size = 16
	point v;
	typedef point non_color_class;

	vert_color() {}
 	vert_color(point const &v_, color_wrapper const &cw) : color_wrapper(cw), v(v_) {}
	vert_color(point const &v_, colorRGBA const &c_)     : v(v_) {set_c4(c_);}
	vert_color(point const &v_, unsigned char const *c_) : v(v_) {c[0]=c_[0]; c[1]=c_[1]; c[2]=c_[2]; c[3]=c_[3];}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void unset_attrs() {}
};


struct vert_norm_color : public vert_norm, public color_wrapper { // size = 28
	typedef vert_norm non_color_class;
	vert_norm_color() {}
	vert_norm_color(vert_norm const &vn, color_wrapper const &cw) : vert_norm(vn), color_wrapper(cw) {}
	vert_norm_color(vert_norm const &vn, colorRGBA const &c_) : vert_norm(vn) {set_c4(c_);}
	vert_norm_color(point const &v_, vector3d const &n_, colorRGBA const     &c_) : vert_norm(v_, n_) {set_c4(c_);}
	vert_norm_color(point const &v_, vector3d const &n_, unsigned char const *c_) : vert_norm(v_, n_) {c[0]=c_[0]; c[1]=c_[1]; c[2]=c_[2]; c[3]=c_[3];}
	void assign(point const &v_, vector3d const &n_, color_wrapper const &cw) {v = v_; n = n_; copy_color(cw);}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_comp_color : public vert_norm_comp, public color_wrapper { // size = 20
	typedef vert_norm_comp non_color_class; // non-compressed type
	vert_norm_comp_color() {}
	vert_norm_comp_color(vert_norm const &vn, color_wrapper const &cw) : vert_norm_comp(vn), color_wrapper(cw) {}
	vert_norm_comp_color(vert_norm_comp const &vn, color_wrapper const &cw) : vert_norm_comp(vn), color_wrapper(cw) {}
	vert_norm_comp_color(point const &v_, vector3d  const &n_, colorRGB  const &c_) : vert_norm_comp(v_, n_) {set_c3(c_);}
	vert_norm_comp_color(point const &v_, vector3d  const &n_, colorRGBA const &c_) : vert_norm_comp(v_, n_) {set_c4(c_);}
	vert_norm_comp_color(point const &v_, vector3d  const &n_, color_wrapper const &cw) : vert_norm_comp(v_, n_), color_wrapper(cw) {}
	vert_norm_comp_color(point const &v_, norm_comp const &n_, color_wrapper const &cw) : vert_norm_comp(v_, n_), color_wrapper(cw) {}
	void assign(point const &v_, vector3d const &n_, unsigned char const *const c_, bool has_alpha=0) {
		v = v_; set_norm(n_); copy_color(c_, has_alpha);
	}
	void assign(point const &v_, char const *const n_, unsigned char const *const c_, bool has_alpha=0) {
		v = v_; set_norm(n_); copy_color(c_, has_alpha);
	}
	void assign(point const &v_, norm_comp const &n_, color_wrapper const &c_) {
		v = v_; set_norm(n_); copy_color(c_);
	}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_tc_color : public vert_norm_tc, public color_wrapper { // size = 36
	typedef vert_norm_tc non_color_class;
	vert_norm_tc_color() {}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, colorRGB const &c_)
		: vert_norm_tc(v_, n_, ts, tt) {set_c3(c_);}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, colorRGBA const &c_)
		: vert_norm_tc(v_, n_, ts, tt) {set_c4(c_);}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, color_wrapper const &cw)
		: vert_norm_tc(v_, n_, ts, tt), color_wrapper(cw) {}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0)
		: vert_norm_tc(v_, n_, ts, tt) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2]; c[3] = (has_alpha ? c_[3] : 255);}
	vert_norm_tc_color(vert_norm const &vn, float ts, float tt, color_wrapper const &cw) : vert_norm_tc(vn, ts, tt), color_wrapper(cw) {}
	vert_norm_tc_color(vert_norm_tc const &vntc, color_wrapper const &cw) : vert_norm_tc(vntc), color_wrapper(cw) {}
	void assign(point const &v_, vector3d const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0) {
		v = v_; n = n_; t[0] = ts; t[1] = tt; copy_color(c_, has_alpha);
	}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_tc_color : public vert_tc_t, public color_wrapper { // size = 24
	typedef vert_tc_t non_color_class;
	vert_tc_color() {}
	vert_tc_color(point const &v_, float ts, float tt, colorRGBA const &c_) : vert_tc_t(v_, ts, tt) {set_c4(c_);}
	vert_tc_color(point const &v_, float ts, float tt, unsigned char const c_[4]) : vert_tc_t(v_, ts, tt) {UNROLL_4X(c[i_] = c_[i_];)}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_comp_tc_color : public vert_norm_comp_tc, public color_wrapper { // size = 28
	typedef vert_norm_tc non_color_class;
	vert_norm_comp_tc_color() {}
	vert_norm_comp_tc_color(vert_norm_comp_tc const &vntc, color_wrapper const &cw) : vert_norm_comp_tc(vntc), color_wrapper(cw) {}
	template<typename T> vert_norm_comp_tc_color(point const &v_, T const &n_, float ts, float tt, color_wrapper const &cw) : vert_norm_comp_tc(v_, n_, ts, tt), color_wrapper(cw) {}
	template<typename T> void assign(point const &v_, T const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0) { // T can be vector3d or norm_comp
		v = v_; set_norm(n_); t[0] = ts; t[1] = tt; copy_color(c_, has_alpha);
	}
	template<typename T> void assign(point const &v_, T const &n_, float ts, float tt, color_wrapper const &cw) { // T can be vector3d or norm_comp
		v = v_; set_norm(n_); t[0] = ts; t[1] = tt; copy_color(cw);
	}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_comp_tc_comp_color : public vert_norm_comp_tc_comp, public color_wrapper { // size = 24
	typedef vert_norm_tc non_color_class;
	vert_norm_comp_tc_comp_color() {}
	vert_norm_comp_tc_comp_color(vert_norm_comp_tc_comp const &vntc, color_wrapper const &cw) : vert_norm_comp_tc_comp(vntc), color_wrapper(cw) {}
	void assign(point const &v_, vector3d const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0) {
		v = v_; set_norm(n_); t[0] = 32767*ts; t[1] = 32767*tt; copy_color(c_, has_alpha);
	}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
};


struct vert_norm_color_tangent : public vert_norm_color {
	vector3d t;

	vert_norm_color_tangent() {}
	vert_norm_color_tangent(point const &v_, vector3d const &n_, colorRGBA const     &c_, vector3d const &t_) : vert_norm_color(v_, n_, c_), t(t_) {}
	vert_norm_color_tangent(point const &v_, vector3d const &n_, unsigned char const *c_, vector3d const &t_) : vert_norm_color(v_, n_, c_), t(t_) {}
};


struct texgen_params_t { // size = 32
	float st[2][4];
	texgen_params_t() {UNROLL_4X(st[0][i_] = st[1][i_] = 0.0;)} // zero initialized
};

// Note: could probably use norm_comp, since used for cubes, cylinder ends, and polygons, but won't save much
struct vert_norm_texp : public vert_norm, public texgen_params_t { // size = 76
	vert_norm_texp() {}
	vert_norm_texp(vert_norm const &vn, texgen_params_t const &tp) : vert_norm(vn), texgen_params_t(tp) {}
	vert_norm_texp(point const &v_, vector3d const &n_, texgen_params_t const &tp) : vert_norm(v_, n_), texgen_params_t(tp) {}
	vert_norm_texp(float x, float y, float z, vector3d const &n_, texgen_params_t const &tp) : vert_norm(point(x, y, z), n_), texgen_params_t(tp) {}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void unset_attrs();
};

void const *ptr_add(void const *p, unsigned off);
void set_vn_ptrs(unsigned stride, bool comp, void const *vbo_ptr_offset);
