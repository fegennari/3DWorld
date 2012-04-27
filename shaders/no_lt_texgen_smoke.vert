uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform mat4 world_space_mvm;

attribute vec4 tex0_s, tex0_t;

varying vec3 eye, vpos, normal, lpos0, vposl; // world space
// epos and eye_norm come from bump_map.vert

void main()
{
	if (use_texgen == 1) {
		setup_texgen(0);
	}
	else if (use_texgen == 2) {
		gl_TexCoord[0].s = dot(gl_Vertex, tex0_s);
		gl_TexCoord[0].t = dot(gl_Vertex, tex0_t);
		gl_TexCoord[0].p = 0.0;
		gl_TexCoord[0].q = 1.0;
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	eye_norm = normalize(gl_NormalMatrix * gl_Normal);
	epos     = gl_ModelViewMatrix * gl_Vertex;

	if (use_world_space_mvm) {
		normal = normalize((transpose(world_space_mvm) * vec4(eye_norm, 1)).xyz);
		mat4 mvm_inv = inverse(world_space_mvm);
		vpos   = (mvm_inv * epos).xyz; // world space
		eye    = mvm_inv[3].xyz; // world space
	}
	else {
		vpos = gl_Vertex.xyz;
		eye  = gl_ModelViewMatrixInverse[3].xyz; // world space
		normal   = normalize(gl_Normal);
	}
	setup_indir_lighting(vpos, normal);

#ifdef USE_BUMP_MAP
	setup_tbn();
#endif

	if (!smoke_enabled) {
		set_fog(); // set standard fog coord
	}
	else if (dynamic_smoke_shadows) {
		lpos0 = (gl_ModelViewMatrixInverse * gl_LightSource[0].position).xyz;
		pt_pair res2 = clip_line(vpos, lpos0, smoke_bb);
		lpos0 = res2.v1;
		vposl = res2.v2;
	}
} 
