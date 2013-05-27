uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform mat4 world_space_mvm;
uniform float tex_scale_s  = 1.0;
uniform float tex_scale_t  = 1.0;
uniform vec3 world_space_offset = vec3(0.0);

attribute vec4 tex0_s, tex0_t;

varying vec3 eye, vpos, normal, lpos0, vposl; // world space
varying vec2 tex_coord; // FIXME: why doesn't gl_TexCoord[0] work?
// epos and eye_norm come from bump_map.vert

void main()
{
	if (use_texgen == 1) {
		setup_texgen0();
	}
	else if (use_texgen == 2) {
		gl_TexCoord[0].s = dot(gl_Vertex, tex0_s);
		gl_TexCoord[0].t = dot(gl_Vertex, tex0_t);
	}
	else if (use_texgen == 3) {
		set_tc0_from_vert_id();
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
		gl_TexCoord[0].st *= vec2(tex_scale_s, tex_scale_t);
	}
	tex_coord     = gl_TexCoord[0].st;
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	eye_norm      = normalize(gl_NormalMatrix * gl_Normal);
	epos          = gl_ModelViewMatrix * gl_Vertex;

	if (use_world_space_mvm) {
		normal = normalize((transpose(world_space_mvm) * vec4(eye_norm, 1)).xyz);
		mat4 mvm_inv = inverse(world_space_mvm);
		vpos   = (mvm_inv * epos).xyz; // world space
		eye    = mvm_inv[3].xyz; // world space
	}
	else {
		vpos   = gl_Vertex.xyz + world_space_offset;
		eye    = gl_ModelViewMatrixInverse[3].xyz; // world space
		normal = normalize(gl_Normal);
	}
	setup_indir_lighting(vpos, normal);
#ifdef USE_BUMP_MAP
	setup_tbn();
#endif

#ifndef SMOKE_ENABLED
#endif
#ifdef DYNAMIC_SMOKE_SHADOWS
	lpos0 = (gl_ModelViewMatrixInverse * gl_LightSource[0].position).xyz;
	pt_pair res2 = clip_line(vpos, lpos0, smoke_bb);
	lpos0 = res2.v1;
	vposl = res2.v2;
#endif
} 
