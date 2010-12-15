uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
varying vec3 eye, vpos, spos, dlpos, normal;
varying vec3 eye_norm;

void main()
{
	if (use_texgen) {
		setup_texgen(0);
	}
	else {
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}	
	gl_Position = ftransform();
	gl_FrontColor = gl_Color;
	normal   = normalize(gl_Normal);
	eye_norm = normalize(gl_NormalMatrix * gl_Normal);
	
	dlpos = gl_Vertex.xyz;
	spos  = gl_Vertex.xyz + (0.25*step_delta)*normal; // move slightly away from the vertex
	eye   = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	
	if (!smoke_enabled) { // set t zero length vector
		vpos = eye; // Note: eye is used for dynamic lights, but vpos is not
		return;
	}
	pt_pair res = clip_line(gl_Vertex.xyz, eye, smoke_bb);
	eye  = res.v1;
	vpos = res.v2;
} 
