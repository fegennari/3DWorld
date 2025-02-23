// 3D World - gl_ext/gl_arb wrapper function declarations
// by Frank Gennari
// 7/31/06
#pragma once

unsigned const PRIMITIVE_RESTART_IX = 0xFFFFFFFF;
unsigned const NUM_TEX_MS_SAMPLES   = 4; // or 8


GLenum get_internal_texture_format(int ncolors, bool compressed=0, bool linear_space=0); // Note: ncolors=2 is unused

inline GLenum get_texture_format(int ncolors) {
	assert(ncolors >= 1 && ncolors <= 4);
	GLenum const formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
	return formats[ncolors-1];
}


void init_glew();

// multitexture prototypes
void bind_texture_tu(unsigned tid, unsigned tu_id);
void bind_texture_tu_def_white_tex(unsigned tid, unsigned tu_id);

// 3D texture prototypes
void bind_3d_texture(unsigned tid);
void setup_3d_texture(unsigned &tid, int filter, int wrap);
unsigned create_3d_texture(unsigned xsz, unsigned ysz, unsigned zsz, unsigned ncomp, vector<unsigned char> const &data, int filter, int wrap, bool compress=0, unsigned bytes_per_pixel=1);
void update_3d_texture(unsigned tid, unsigned xoff, unsigned yoff, unsigned zoff, unsigned xsz, unsigned ysz, unsigned zsz,
					   unsigned ncomp, unsigned char const *const data);

// gl_ext_arb
unsigned create_vbo();
void bind_vbo(unsigned vbo, bool is_index=0);
void delete_vbo(unsigned vbo);
void upload_vbo_data(void const *const data, size_t size, bool is_index=0, int dynamic_level=0);
void upload_vbo_sub_data(void const *const data, int offset, size_t size, bool is_index=0);
void upload_vbo_sub_data_no_sync(void const *data, unsigned start_byte, unsigned size_bytes, bool is_index=0);
void bind_ubo(unsigned ubo);
void bind_ubo_base(unsigned ubo, unsigned index);
void upload_ubo_data(void const *const data, size_t size, int dynamic_level=0);
void upload_ubo_sub_data(void const *const data, int offset, size_t size);
unsigned create_vao();
void bind_vao(unsigned vao);
void delete_vao(unsigned vao);
void create_fbo(unsigned &fbo_id, unsigned depth_tid, bool is_depth_fbo=0, bool multisample=0, bool is_array=0, unsigned *layer=nullptr);
void enable_fbo(unsigned &fbo_id, unsigned tid,       bool is_depth_fbo=0, bool multisample=0, bool is_array=0, unsigned *layer=nullptr);
void bind_fbo(unsigned fbo_id);
void disable_fbo();
void free_fbo(unsigned &fbo_id);
void bind_pbo(unsigned pbo_id);
unsigned create_depth_render_buffer(unsigned xsize, unsigned ysize, bool multisample=0);
void disable_and_free_render_buffer(unsigned &render_buffer);
bool gen_mipmaps(unsigned dim=2);
void enable_instancing_for_shader_loc(int loc);
void disable_instancing_for_shader_loc(int loc);

inline void delete_and_zero_vbo(unsigned &vbo) {delete_vbo(vbo); vbo = 0;}
inline void delete_and_zero_vao(unsigned &vao) {delete_vao(vao); vao = 0;}
inline void check_bind_vbo(unsigned vbo, bool is_index=0) {assert(vbo); bind_vbo(vbo, is_index);}
inline void check_bind_ubo(unsigned ubo) {assert(ubo); bind_ubo(ubo);}
inline void check_bind_vao(unsigned vao) {assert(vao); bind_vao(vao);}

void setup_stencil_buffer_write();
void end_stencil_write();


// templated vbo management utility functions/classes

template<typename T> void upload_to_vbo(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0, int dynamic_level=0) {
	check_bind_vbo(vbo, is_index);
	upload_vbo_data(data.data(), data.size()*sizeof(T), is_index, dynamic_level);
	if (end_with_bind0) {bind_vbo(0, is_index);}
}
template<typename T> void upload_vector_to_vbo(vector<T> const &data, bool is_index=0, int offset=0) {
	upload_vbo_sub_data(data.data(), offset, data.size()*sizeof(T), is_index);
}

template<typename T> void upload_to_ubo(unsigned &vbo, vector<T> const &data, bool end_with_bind0=0, int dynamic_level=0) {
	check_bind_ubo(vbo);
	upload_ubo_data(data.data(), data.size()*sizeof(T), dynamic_level);
	if (end_with_bind0) {bind_ubo(0);}
}
template<typename T> void upload_vector_to_ubo(vector<T> const &data) {
	upload_ubo_sub_data(data.data(), 0, data.size()*sizeof(T)); // offset=0
}

template<typename T> bool create_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0, bool end_with_bind0=0, int dynamic_level=0) {
	if (vbo) return 0; // already uploaded
	vbo = create_vbo();
	upload_to_vbo(vbo, data, is_index, end_with_bind0, dynamic_level);
	return 1;
}
template<typename T> bool create_ubo_and_upload(unsigned &ubo, vector<T> const &data, bool end_with_bind0=0, int dynamic_level=0) {
	if (ubo) return 0; // already uploaded
	ubo = create_vbo(); // same as VBO
	upload_to_ubo(ubo, data, end_with_bind0, dynamic_level);
	return 1;
}

template<typename T> void create_bind_vbo_and_upload(unsigned &vbo, vector<T> const &data, bool is_index=0, int dynamic_level=0) {
	if (!create_vbo_and_upload(vbo, data, is_index, 0, dynamic_level)) {bind_vbo(vbo, is_index);}
}

inline void create_vbo_with_null_data(unsigned &vbo, size_t size, bool is_index=0, int dynamic_level=0) {
	vbo = create_vbo();
	bind_vbo(vbo, is_index);
	upload_vbo_data(NULL, size, is_index, dynamic_level);
}

template<typename T> void update_indices(unsigned ivbo, vector<T> const &indices, bool end_with_bind_0=1) {
	check_bind_vbo(ivbo, 1);
	upload_vector_to_vbo(indices, 1);
	if (end_with_bind_0) {bind_vbo(0, 1);}
}


struct vbo_wrap_t { // Note: not for use with index vbo

	unsigned vbo=0;

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

struct ubo_wrap_t { // uniform buffer object

	unsigned ubo=0;

	bool ubo_valid() const {return (ubo > 0);}
	void clear() {delete_and_zero_vbo(ubo);} // same as VBO
	void allocate_with_size(unsigned size, int dynamic_level=0);
	template<typename vert_type_t>
	void create_and_upload(vector<vert_type_t> const &data, int dynamic_level=0, bool end_with_bind0=0) {
		if (!ubo) {create_ubo_and_upload(ubo, data, end_with_bind0, dynamic_level);}
	}
	void bind_base(unsigned index) {bind_ubo_base(ubo, index);}
	void pre_render() const {check_bind_ubo(ubo);}
	static void post_render() {bind_ubo(0);}
};

struct vao_wrap_t {

	unsigned vao=0;

	bool is_valid() const {return (vao != 0);}
	void clear() {delete_and_zero_vao(vao);}

	void ensure_vao_bound() {
		if (!vao) {vao = create_vao();}
		enable_vao();
	}
	void enable_vao() const {check_bind_vao(vao);}
	static void disable_vao() {bind_vao(0);}

	template<typename vert_type_t> void create_from_vbo(vbo_wrap_t const &vbo, bool setup_pointers=0, bool always_bind=0) {
		if (vao) {if (always_bind) {enable_vao();} return;} // already set
		ensure_vao_bound();
		vbo.pre_render();
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
	}
};

struct indexed_vbo_manager_t : public vbo_wrap_t {

	unsigned ivbo=0, gpu_mem=0; // aka EBO

	bool ivbo_valid() const {return (ivbo > 0);}
	void reset_vbos_to_zero() {vbo = ivbo = gpu_mem = 0;}

	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata, int dynamic_level=0, bool end_with_bind0=0) {
		vbo_wrap_t::create_and_upload(data, dynamic_level, end_with_bind0);
		if (!vbo ) {gpu_mem += data.size() *sizeof(vert_type_t );}
		if (!ivbo && !idata.empty()) {create_vbo_and_upload(ivbo, idata, 1, end_with_bind0, dynamic_level); gpu_mem += idata.size()*sizeof(index_type_t);}
	}
	template<typename vert_type_t, typename index_type_t>
	void create_and_upload_vbo_with_vao(vector<vert_type_t> const &data, vector<index_type_t> const &idata, vao_wrap_t &vao_wrap, int dynamic_level=0, bool setup_pointers=0) {
		if (vao_wrap.vao) return; // already set
		vao_wrap.ensure_vao_bound();
		if (vbo) {indexed_vbo_manager_t::pre_render(!idata.empty());} else {indexed_vbo_manager_t::create_and_upload(data, idata, dynamic_level);}
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
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

struct vao_manager_t : public vbo_wrap_t, public vao_wrap_t {

	template<typename vert_type_t> void create_and_upload(vector<vert_type_t> const &data, int dynamic_level=0, bool setup_pointers=0) {
		if (vao) {assert(vbo); return;} // already set
		ensure_vao_bound();
		if (vbo) {vbo_wrap_t::pre_render();} else {vbo_wrap_t::create_and_upload(data, dynamic_level);}
		if (setup_pointers) {vert_type_t::set_vbo_arrays();}
	}
	void clear() {vbo_wrap_t::clear(); vao_wrap_t::clear();}
	void pre_render(bool do_bind_vbo=0) const {enable_vao(); if (do_bind_vbo) {vbo_wrap_t::pre_render();}}
	static void post_render() {disable_vao(); vbo_wrap_t::post_render();}
};

struct indexed_vao_manager_t : public indexed_vbo_manager_t, public vao_wrap_t {

	void reset_vbos_to_zero() {indexed_vbo_manager_t::reset_vbos_to_zero(); vao = 0;} // virtual?
	void clear_vbos() {indexed_vbo_manager_t::clear_vbos(); vao_wrap_t::clear();} // and VAO

	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata, int dynamic_level=0, bool setup_pointers=0) {
		create_and_upload_vbo_with_vao(data, idata, *this, dynamic_level, setup_pointers);
	}
	void pre_render(bool using_index=1, bool do_bind_vbo=0) const {
		enable_vao();
		if (do_bind_vbo) {indexed_vbo_manager_t::pre_render(using_index);}
	}
	static void post_render() {disable_vao(); indexed_vbo_manager_t::post_render();}
};

struct indexed_vao_manager_with_shadow_t : public indexed_vbo_manager_t {
protected:
	vao_wrap_t vaos[2]; // {regular, shadow}
public:
	bool is_vao_setup(bool shadow) const {return vaos[shadow].is_valid();}
	void reset_vbos_to_zero() {indexed_vbo_manager_t::reset_vbos_to_zero(); vaos[0].vao = vaos[1].vao = 0;} // virtual?
	void clear_vaos() {vaos[0].clear(); vaos[1].clear();}
	void clear_vbos() {indexed_vbo_manager_t::clear_vbos(); clear_vaos();} // and VAOs

	template<typename vert_type_t, typename index_type_t>
	void create_and_upload(vector<vert_type_t> const &data, vector<index_type_t> const &idata, bool shadow, int dynamic_level=0, bool setup_pointers=0) {
		create_and_upload_vbo_with_vao(data, idata, vaos[shadow], dynamic_level, setup_pointers);
	}
	template<typename vert_type_t> void create_from_vbo(bool shadow, bool setup_pointers=0, bool always_bind=0) {
		vaos[shadow].create_from_vbo<vert_type_t>(*this, setup_pointers, always_bind);
	}
	void pre_render(bool shadow, bool using_index=1, bool do_bind_vbo=0) const {
		vaos[shadow].enable_vao();
		if (do_bind_vbo) {indexed_vbo_manager_t::pre_render(using_index);}
	}
	static void post_render() {bind_vao(0); indexed_vbo_manager_t::post_render();}
};

template<unsigned N> struct indexed_vao_multi_manager_t : public indexed_vbo_manager_t { // unused

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
		create_and_upload_vbo_with_vao(data, idata, vaos[ix], dynamic_level, setup_pointers);
	}
	void enable_vao(unsigned ix) const {assert(ix < N); vaos[ix].enable_vao();}
	void pre_render(unsigned ix, bool using_index=1) const {enable_vao(ix); indexed_vbo_manager_t::pre_render(using_index);}
};


class subdiv_sphere_drawer_t : public indexed_vao_manager_t {
protected:
	unsigned nverts=0, nindices=0;
};
struct icosphere_drawer_t : public subdiv_sphere_drawer_t {
	icosphere_drawer_t(unsigned ndiv);
	void draw() const;
};

template<typename T> class subdiv_sphere_manager_t {
	typedef map<unsigned, T> ndiv_sphere_map_t;
	ndiv_sphere_map_t cached;
public:
	void draw_sphere(unsigned ndiv);
	void clear();
};

typedef subdiv_sphere_manager_t<icosphere_drawer_t> icosphere_manager_t;


class vbo_ring_buffer_t : public vbo_wrap_t {

	unsigned init_size, size, pos=0;
	bool is_index;

	void ensure_vbo(unsigned min_size);
public:
	vbo_ring_buffer_t(unsigned init_size_, bool is_index_=0) : init_size(init_size_), size(init_size), is_index(is_index_) {}
	void clear() {vbo_wrap_t::clear(); size = init_size; pos = 0;}
	unsigned get_alloced_size() const {return (vbo_valid() ? size : 0);}

	template<typename T> void *add_verts_bind_vbo(vector<T> const &v) {
		assert(!v.empty());
		return add_verts_bind_vbo(&v.front(), v.size()*sizeof(T));
	}
	bool has_space_for(unsigned size_bytes) const {return (pos + size_bytes <= size);}
	void mark_as_filled() {pos = size;}
	void orphan_buffer_and_reset();
	void const *add_verts_bind_vbo(void const *const v, unsigned size_bytes);
};

inline void align_vbo_ptr(unsigned &pos) {if (pos & 15) {pos = (pos + 16) & (~15);}}


void const *get_dynamic_vbo_ptr(void const *const verts, unsigned size_bytes);
void bind_dynamic_vbo();
void ensure_texture_loaded(unsigned &tid, unsigned txsize, unsigned tysize, bool mipmap, bool nearest, bool multisample=0);
void build_texture_mipmaps(unsigned tid, unsigned dim);


struct texture_pair_t {

	unsigned tids[2]={}; // color, normal
	bool multisample;

	texture_pair_t(bool multisample_=0) : multisample(multisample_) {}
	bool is_valid() const {return (tids[0] > 0 && tids[1] > 0);}
	void free_context();
	void bind_texture() const;
	void ensure_tid(unsigned tsize, bool mipmap);
	bool operator==(texture_pair_t const &tp) const {return (tids[0] == tp.tids[0] && tids[1] == tp.tids[1]);}
	bool operator!=(texture_pair_t const &tp) const {return !operator==(tp);}
	bool operator< (texture_pair_t const &tp) const {return ((tids[0] == tp.tids[0]) ? (tids[1] < tp.tids[1]) : (tids[0] < tp.tids[0]));}
};

struct texture_atlas_t { // unused

	unsigned tid=0, nx=1, ny=1;
	bool multisample;

	texture_atlas_t(bool multisample_=0) : multisample(multisample_) {}
	texture_atlas_t(unsigned nx_, unsigned ny_, bool multisample_=0) : nx(nx_), ny(ny_), multisample(multisample_) {}
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
	virtual ~render_to_texture_t() {}
	void render(texture_pair_t &tpair, float xsize, float ysize, point const &center, vector3d const &view_dir,
		colorRGBA const &bkg_color, bool use_depth_buffer, bool mipmap);
	virtual void draw_geom(bool is_normal_pass) = 0;
};

void set_temp_clear_color(colorRGBA const &clear_color, bool clear_depth=0, bool clear_stencil=0);

template< typename T > void upload_to_dynamic_vbo(vector<T> const &v) {
	T::set_vbo_arrays(1, get_dynamic_vbo_ptr(&v.front(), v.size()*sizeof(T)));
}

void bind_ssbo(unsigned ssbo=0);
unsigned create_ssbo(unsigned data_sz=0, void const *const data=nullptr);
void ensure_ssbo(unsigned &ssbo, unsigned data_sz, void const *const data);
void update_ssbo(unsigned  ssbo, unsigned data_sz, void const *const data);



class query_perf_timer_t {
	unsigned query_id=0;
public:
	query_perf_timer_t();
	~query_perf_timer_t();
	void time_query();
	GLuint64 get_query();
};

class gpu_timer_t {
	query_perf_timer_t time_query;
	GLint64 timer=0;
public:
	gpu_timer_t();
	void show();
};

GLint64 get_timestamp();

struct DrawElementsIndirectCommand { // used with glMultiDrawElementsIndirect()
	uint32_t count=0, instanceCount=0, firstIndex=0;
	int32_t  baseVertex=0;
	uint32_t baseInstance=0; // OpenGL >= 4.2
};

inline int get_2d_texture_target(bool is_array=0, bool multisample=0) {
	return (is_array ? (multisample ? GL_TEXTURE_2D_MULTISAMPLE_ARRAY : GL_TEXTURE_2D_ARRAY) : (multisample ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D));
}

class DebugScope {
	static GLuint global_scope_depth;
	const GLuint scope_depth;
public:
	DebugScope(std::string const &name) : scope_depth(global_scope_depth++) {
		glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, scope_depth, name.size(), name.data());
	}
	~DebugScope() {
		glPopDebugGroup();
		global_scope_depth--;
	}
};


