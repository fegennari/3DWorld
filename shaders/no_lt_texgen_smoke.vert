uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;

attribute float shadow_val; // sending as int doesn't work?

varying vec3 eye, vpos, spos, dlpos, normal, lpos0, vposl; // world space
varying vec3 eye_norm;
varying float light_scale[8];

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position   = ftransform();
	gl_FrontColor = gl_Color;
	normal   = normalize(gl_Normal);
	eye_norm = normalize(gl_NormalMatrix * gl_Normal);

	dlpos    = gl_Vertex.xyz;
	spos     = gl_Vertex.xyz + (0.25*step_delta)*normal; // move slightly away from the vertex
	eye      = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0  )).xyz; // world space
	int shadow_bits = int(round(shadow_val));

	for (uint i = 0; i < 8; ++i) {
		light_scale[i] = (((shadow_bits & (1 << i)) == 0) ? 1.0 : 0.0);
	}
	if (!smoke_enabled) { // set t zero length vector
		vpos = eye; // Note: eye is used for dynamic lights, but vpos is not
		set_fog(); // set standard fog coord
		return;
	}
	pt_pair res = clip_line(gl_Vertex.xyz, eye, smoke_bb);
	eye  = res.v1;
	vpos = res.v2;

	if (dynamic_smoke_shadows) {
		lpos0 = (gl_ModelViewMatrixInverse * gl_LightSource[0].position).xyz;
		pt_pair res2 = clip_line(gl_Vertex.xyz, lpos0, smoke_bb);
		lpos0 = res2.v1;
		vposl = res2.v2;
	}
} 
