uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform mat4 world_space_mvm;
uniform float tex_scale_s  = 1.0;
uniform float tex_scale_t  = 1.0;
uniform vec3 world_space_offset = vec3(0.0);
uniform vec3 sun_pos; // used for dynamic smoke shadows line clipping

attribute vec4 tex0_s, tex0_t;

varying vec3 vpos, normal, lpos0, vposl; // world space
// epos, eye_norm, and tex_coord come from bump_map.vert

void main()
{
	if (use_texgen == 1) {
		setup_texgen_st();
		tex_coord = tc;
	}
	else if (use_texgen == 2) {
		tex_coord = vec2(dot(fg_Vertex, tex0_s), dot(fg_Vertex, tex0_t));
	}
	else if (use_texgen == 3) {
		set_tc0_from_vert_id();
		tex_coord = tc;
	}
	else {
		tex_coord = fg_TexCoord * vec2(tex_scale_s, tex_scale_t);
	}
	gl_FrontColor = fg_Color;
	epos          = gl_ModelViewMatrix * fg_Vertex;
	gl_Position   = gl_ProjectionMatrix * epos;

	if (use_world_space_mvm) {
		eye_norm = normalize(gl_NormalMatrix * fg_Normal);
		normal   = normalize((transpose(world_space_mvm) * vec4(eye_norm, 1)).xyz);
		vpos     = (inverse(world_space_mvm) * epos).xyz; // world space
	}
	else {
		eye_norm = normalize(mat3(gl_ModelViewMatrix) * fg_Normal); // Note: avoids the gl_NormalMatrix upload
		vpos     = fg_Vertex.xyz + world_space_offset;
		normal   = normalize(fg_Normal);
	}
#ifdef USE_BUMP_MAP
	setup_tbn();
#endif

#ifdef DYNAMIC_SMOKE_SHADOWS
	pt_pair res2 = clip_line(vpos, sun_pos, smoke_bb);
	lpos0 = res2.v1;
	vposl = res2.v2;
#endif
} 
