uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float half_dxy;
uniform float indir_vert_offset = 0.25;

attribute vec4 tex0_s, tex0_t;

varying vec3 eye, vpos, spos, normal, lpos0, vposl; // world space
varying vec3 eye_norm;
varying vec4 epos;

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
	normal   = normalize(gl_Normal);
	eye_norm = normalize(gl_NormalMatrix * gl_Normal);
	epos     = gl_ModelViewMatrix * gl_Vertex;
	vpos     = gl_Vertex.xyz;
	spos     = gl_Vertex.xyz + (indir_vert_offset*half_dxy)*normal; // move slightly away from the vertex
	eye      = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space

#ifdef USE_BUMP_MAP
	setup_tbn(eye_norm, epos.xyz);
#endif

	if (!smoke_enabled) {
		set_fog(); // set standard fog coord
		return;
	}
	if (dynamic_smoke_shadows) {
		lpos0 = (gl_ModelViewMatrixInverse * gl_LightSource[0].position).xyz;
		pt_pair res2 = clip_line(gl_Vertex.xyz, lpos0, smoke_bb);
		lpos0 = res2.v1;
		vposl = res2.v2;
	}
} 
