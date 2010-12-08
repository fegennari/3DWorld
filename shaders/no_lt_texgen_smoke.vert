uniform float smoke_bb[6]; // x1,x2,y1,y2,z1,z2
uniform float step_delta;
varying vec3 eye, vpos, spos, dlpos, normal;

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
	normal = normalize(gl_Normal);
	
	dlpos = gl_Vertex.xyz;
	spos = gl_Vertex.xyz + (0.25*step_delta)*normal; // move slightly away from the vertex
	
	if (!smoke_enabled) {
		eye  = vec3(0,0,0);
		vpos = vec3(0,0,0);
		return;
	}
	vec3 v2 = (gl_ModelViewMatrixInverse * vec4(0.0, 0.0, 0.0, 1.0)).xyz; // world space
	pt_pair res = clip_line(gl_Vertex.xyz, v2, smoke_bb);
	eye  = res.v1;
	vpos = res.v2;
} 
