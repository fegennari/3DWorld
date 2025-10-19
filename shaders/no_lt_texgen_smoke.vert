uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform mat4 fg_ViewMatrix;
uniform float tex_scale_s   = 1.0;
uniform float tex_scale_t   = 1.0;
uniform float tc_texgen_mix = 0.0;
uniform vec4 world_space_offset   = vec4(0.0); // {x, y, z, rot_angle}
uniform float vertex_offset_scale = 0.0; // hack to make vertex_offset ignored when unused/unset

in vec4 tex0_s, tex0_t;
in vec3 vertex_offset; // not always used

out vec3 vpos, normal; // world space
out vec4 epos;
out vec3 eye_norm;
// tc comes from texture_gen.part.vert

#ifdef ENABLE_CLIP_PLANE
uniform vec4 clip_plane = vec4(0.0);
out float gl_ClipDistance[1];
#endif

void main() {
	if      (use_texgen == 1) {setup_texgen_st();}
	else if (use_texgen == 2) {tc = vec2(dot(fg_Vertex, tex0_s), dot(fg_Vertex, tex0_t));}
	else if (use_texgen == 3) {set_tc0_from_vert_id();}
	else if (use_texgen == 4) {set_bent_quad_tc0_from_vert_id();}
	else if (use_texgen == 5) {setup_texgen_st(); tc = mix(tc, fg_TexCoord, tc_texgen_mix);}
	else if (use_texgen == 6) {setup_texgen_st_no_xy_cancel(); tc = mix(tc, fg_TexCoord, tc_texgen_mix);}
	else                      {tc = fg_TexCoord * vec2(tex_scale_s, tex_scale_t);}

	vec4 vertex    = vec4((vertex_offset_scale*vertex_offset), 0.0) + fg_Vertex;
	vec3 normal_in = fg_Normal;
#ifdef ENABLE_VERTEX_ANIMATION
	apply_vertex_animation(vertex, normal_in, tc);
#endif
	add_leaf_wind(vertex);
	epos        = fg_ModelViewMatrix * vertex;
	gl_Position = fg_ProjectionMatrix * epos;
	fg_Color_vf = fg_Color;
	eye_norm    = normalize(fg_NormalMatrix * normal_in);

	if (use_fg_ViewMatrix) {
		normal = normalize((transpose(fg_ViewMatrix) * vec4(eye_norm, 1)).xyz);
		vpos   = (inverse(fg_ViewMatrix) * epos).xyz; // world space
	}
	else {
		vpos   = vertex.xyz + world_space_offset.xyz; // Note: rotation not supported here
		normal = normalize(normal_in);
	}
#ifdef USE_BUMP_MAP
	setup_tbn();
#endif

#ifdef ENABLE_CLIP_PLANE
	gl_ClipDistance[0] = dot(vec4(vpos, 1.0), clip_plane);
#endif
} 
