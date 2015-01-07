// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06

#ifndef _GL_EXT_ARB_H_
#define _GL_EXT_ARB_H_

unsigned const PRIMITIVE_RESTART_IX = 0xFFFFFFFF;


inline GLenum get_internal_texture_format(int ncolors, bool compressed=0) { // Note: ncolors=2 is unused
	GLenum const cformats[4] = {GL_COMPRESSED_RED, GL_COMPRESSED_RG, GL_COMPRESSED_RGB, GL_COMPRESSED_RGBA};
	GLenum const formats [4] = {GL_R8, GL_RG8, GL_RGB8, GL_RGBA8};
	return (compressed ? cformats : formats)[ncolors-1];
}

inline GLenum get_texture_format(int ncolors) {
	assert(ncolors >= 1 && ncolors <= 4);
	GLenum const formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
	return formats[ncolors-1];
}


void init_glew();

// multitexture prototypes
void set_active_texture(unsigned tu_id);
void select_multitex(int id, unsigned tu_id, bool reset=1);

// 3D texture prototypes
void bind_3d_texture(unsigned tid);
unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data, int filter, int wrap, bool compress=0);
void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data);
void set_3d_texture_as_current(unsigned tid, unsigned tu_id);

// gl_ext_arb
unsigned create_vbo();
void bind_vbo(unsigned vbo, bool is_index=0);
void delete_vbo(unsigned vbo);
void upload_vbo_data(void const *const data, size_t size, bool is_index=0, int dynamic_level=0);
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index=0);
void upload_vbo_sub_data_no_sync(void const *data, unsigned start_byte, unsigned size_bytes, bool is_index=0);
unsigned create_vao();
void bind_vao(unsigned vao);
void delete_vao(unsigned vao);
void create_fbo(unsigned &fbo_id, unsigned depth_tid, bool is_depth_fbo);
void enable_fbo(unsigned &fbo_id, unsigned tid, bool is_depth_fbo);
void disable_fbo();
void free_fbo(unsigned &fbo_id);
unsigned create_depth_render_buffer(unsigned xsize, unsigned ysize);
void disable_and_free_render_buffer(unsigned &render_buffer);
bool gen_mipmaps(unsigned dim=2);

inline void delete_and_zero_vbo(unsigned &vbo) {delete_vbo(vbo); vbo = 0;}
inline void delete_and_zero_vao(unsigned &vao) {delete_vao(vao); vao = 0;}
inline void check_bind_vbo(unsigned vbo, bool is_index=0) {assert(vbo); bind_vbo(vbo, is_index);}
inline void check_bind_vao(unsigned vao) {assert(vao); bind_vao(vao);}


// templated vbo management utility functions/classes

template<typename T> void upload_to_vbo(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0, int dynamic_level=0) {
	check_bind_vbo(vbo, is_index);
	upload_vbo_data(&data.front(), data.size()*sizeof(T), is_index, dynamic_level);
	if (end_with_bind0) {bind_vbo(0, is_index);}
}

template<typename T> bool create_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0, int dynamic_level=0) {
	if (vbo) return 0; // already uploaded
	vbo = create_vbo();
	upload_to_vbo(vbo, data, is_index, end_with_bind0, dynamic_level);
	return 1;
}

template<typename T> void create_bind_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0, int dynamic_level=0) {
	if (!create_vbo_and_upload(vbo, data, is_index, 0, dynamic_level)) {bind_vbo(vbo, is_index);}
}

inline void create_vbo_with_null_data(unsigned &vbo, size_t size, bool is_index=0, int dynamic_level=0) {
	vbo = create_vbo();
	bind_vbo(vbo, is_index);
	upload_vbo_data(NULL, size, is_index);
}


struct vbo_wrap_t { // Note: not for use with index vbo

	unsigned vbo;

	vbo_wrap_t() : vbo(0) {}
	bool vbo_valid() const {return (vbo > 0);}
	void clear() {delete_and_zero_vbo(vbo);}
	void clear_vbo() {clear();} // alias for clear()
	template<typename vert_type_t>
	void create_and_upload(vector<vert_type_t> const &data, int dynamic_level=0, bool end_with_bind0=0) {
		if (!vbo) {create_vbo_and_upload(vbo, data, 0, end_with_bind0, dynamic_level);}
	}
	void pre_render() const {check_bind_vbo(vbo);}
	static void post_render() {bind_vbo(0);}
};


struct indexed_vbo_manager_t : public vbo_wrap_t {

	unsigned ivbo, gpu_mem;

	indexed_vbo_manager_t() : ivbo(0), gpu_mem(0) {}
	void reset_vbos_to_zero() {vbo = ivbo = gpu_mem = 0;}

	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata, int dynamic_level=0, bool end_with_bind0=0) {
		vbo_wrap_t::create_and_upload(data, dynamic_level, end_with_bind0);
		if (!vbo ) {gpu_mem += data.size() *sizeof(vert_type_t );}
		if (!ivbo) {create_vbo_and_upload(ivbo, idata, 1, end_with_bind0, dynamic_level); gpu_mem += idata.size()*sizeof(index_type_t);}
	}
	void clear_vbos() {
		vbo_wrap_t::clear();
		delete_vbo(ivbo);
		reset_vbos_to_zero();
	}
	void pre_render(bool using_index=1) const {
		assert(ivbo || !using_index);
		vbo_wrap_t::pre_render();
		bind_vbo(ivbo, 1);
	}
	static void post_render() {
		vbo_wrap_t::post_render();
		bind_vbo(0, 1);
	}
};


struct vao_wrap_t {

	unsigned vao;

	vao_wrap_t() : vao(0) {}
	void clear() {delete_and_zero_vao(vao);}
	
	void ensure_vao_bound() {
		if (!vao) {vao = create_vao();}
		enable_vao();
	}
	void enable_vao() const {check_bind_vao(vao);}
	static void disable_vao() {bind_vao(0);}
};


struct vao_manager_t : public vbo_wrap_t, public vao_wrap_t {

	template<typename vert_type_t>
	void create_and_upload(vector<vert_type_t> const &data, int dynamic_level=0, bool setup_pointers=0) {
		if (vao) {assert(vbo); return;} // already set
		ensure_vao_bound();
		vbo_wrap_t::create_and_upload(data, dynamic_level);
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
	}
	void clear() {vbo_wrap_t::clear(); vao_wrap_t::clear();}
	void pre_render() const {enable_vao(); vbo_wrap_t::pre_render();}
	static void post_render() {disable_vao(); indexed_vbo_manager_t::post_render();}
};


struct indexed_vao_manager_t : public indexed_vbo_manager_t, public vao_wrap_t {

	void reset_vbos_to_zero() {indexed_vbo_manager_t::reset_vbos_to_zero(); vao = 0;} // virtual?
	void clear_vbos() {indexed_vbo_manager_t::clear_vbos(); vao_wrap_t::clear();}

	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata, int dynamic_level=0, bool setup_pointers=0) {
		if (vao) return; // already set
		ensure_vao_bound();
		indexed_vbo_manager_t::create_and_upload(data, idata, dynamic_level);
		//indexed_vbo_manager_t::pre_render(1); // binds vbo/ivbo - may not be necessary
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
	}
	void pre_render(bool using_index=1) const {enable_vao(); indexed_vbo_manager_t::pre_render(using_index);}
	static void post_render() {disable_vao(); indexed_vbo_manager_t::post_render();}
};


template<unsigned N> struct indexed_vao_multi_manager_t : public indexed_vbo_manager_t {

	vao_wrap_t vaos[N];

	void clear_vbos() {
		indexed_vbo_manager_t::clear_vbos();
		for (unsigned i = 0; i < N; ++i) {vaos[i].clear();}
	}
	void ensure_vao_bound(unsigned ix) {
		assert(ix < N);
		vaos[ix].ensure_vao_bound();
	}
	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(unsigned ix, vector<vert_type_t> const &data, vector<index_type_t> const &idata, int dynamic_level=0, bool setup_pointers=0) {
		assert(ix < N);
		if (vao[ix]) return; // already set
		ensure_vao_bound(ix);
		indexed_vbo_manager_t::create_and_upload(data, idata, dynamic_level);
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
	}
	void enable_vao(unsigned ix) const {assert(ix < N); vaos[ix].enable_vao();}
	void pre_render(unsigned ix, bool using_index=1) const {enable_vao(ix); indexed_vbo_manager_t::pre_render(using_index);}
};


class cube_map_sphere_drawer_t : public indexed_vao_manager_t {
	unsigned nverts, nindices;

public:
	cube_map_sphere_drawer_t(unsigned ndiv);
	void draw() const;
};

class cube_map_sphere_manager_t {
	typedef map<unsigned, cube_map_sphere_drawer_t> ndiv_sphere_map_t;
	ndiv_sphere_map_t cached;

public:
	void draw_sphere(unsigned ndiv);
	void clear();
};


class vbo_ring_buffer_t : public vbo_wrap_t {

	unsigned init_size, size, pos;
	bool is_index;

	void ensure_vbo(unsigned min_size);
public:
	vbo_ring_buffer_t(unsigned init_size_, bool is_index_=0) : init_size(init_size_), size(init_size), pos(0), is_index(is_index_) {}
	void clear() {vbo_wrap_t::clear(); size = init_size; pos = 0;}

	template<typename T> void *add_verts_bind_vbo(vector<T> const &v) {
		assert(!v.empty());
		return add_verts_bind_vbo(&v.front(), v.size()*sizeof(T));
	}
	void const *add_verts_bind_vbo(void const *const v, unsigned size_bytes);
};

inline void align_vbo_ptr(unsigned &pos) {if (pos & 15) {pos = (pos + 16) & (~15);}}


void const *get_dynamic_vbo_ptr(void const *const verts, unsigned size_bytes);
void ensure_texture_loaded(unsigned &tid, unsigned txsize, unsigned tysize, bool mipmap, bool nearest);
void build_texture_mipmaps(unsigned tid, unsigned dim);


struct texture_pair_t {

	unsigned tids[2]; // color, normal

	texture_pair_t() {tids[0] = tids[1] = 0;}
	bool is_valid() const {return (tids[0] > 0 && tids[1] > 0);}
	void free_context();
	void bind_texture() const;
	void ensure_tid(unsigned tsize, bool mipmap);
	bool operator==(texture_pair_t const &tp) const {return (tids[0] == tp.tids[0] && tids[1] == tp.tids[1]);}
	bool operator!=(texture_pair_t const &tp) const {return !operator==(tp);}
	bool operator< (texture_pair_t const &tp) const {return ((tids[0] == tp.tids[0]) ? (tids[1] < tp.tids[1]) : (tids[0] < tp.tids[0]));}
};


struct texture_atlas_t {

	unsigned tid, nx, ny;

	texture_atlas_t(unsigned nx_=1, unsigned ny_=1) : tid(0), nx(nx_), ny(ny_) {}
	bool is_valid() const {return (tid > 0);}
	void free_context();
	void bind_texture() const;
	void ensure_tid(unsigned base_tsize, bool mipmap);
	bool operator==(texture_atlas_t const &tp) const {return (tid == tp.tid && nx == tp.nx && ny == tp.ny);} // do we need to compare nx and ny?
	bool operator!=(texture_atlas_t const &tp) const {return !operator==(tp);}
	bool operator< (texture_atlas_t const &tp) const {return (tid < tp.tid);}
};


class render_to_texture_t {

	unsigned tsize;

	void pre_render(float xsize, float ysize, unsigned nx, unsigned ny, point const &center, vector3d const &view_dir) const;
	static void post_render();

public:
	render_to_texture_t(unsigned tsize_) : tsize(tsize_) {}
	virtual ~render_to_texture_t() {free_context();}
	virtual void free_context() {} // nothing to do here
	void render(texture_pair_t &tpair, float xsize, float ysize, point const &center, vector3d const &view_dir,
		colorRGBA const &bkg_color, bool use_depth_buffer, bool mipmap);
	void render(texture_atlas_t &atlas, float xsize, float ysize, point const &center, vector3d const &view_dir,
		colorRGBA const &bkg_color, bool use_depth_buffer, bool mipmap);
	virtual void draw_geom(bool is_normal_pass) = 0;
};

void set_temp_clear_color(colorRGBA const &clear_color);

template< typename T > void upload_to_dynamic_vbo(vector<T> const &v) {
	T::set_vbo_arrays(1, get_dynamic_vbo_ptr(&v.front(), v.size()*sizeof(T)));
}


class query_perf_timer_t {
	unsigned query_id;

public:
	query_perf_timer_t();
	~query_perf_timer_t();
	void time_query();
	GLuint64 get_query();
};

class gpu_timer_t {
	query_perf_timer_t time_query;
	GLint64 timer;

public:
	gpu_timer_t();
	void show();
};

GLint64 get_timestamp();


#endif // _GL_EXT_ARB_H_

